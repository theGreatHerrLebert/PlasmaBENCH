"""Build a TimSim from_findings seed CSV from a real search-engine report.

Reads a DIA-NN or FragPipe report for ONE PYE1 sample and emits the 6-column
findings contract that TimSim's `from_findings` mode consumes:

    protein, sequence, intensity   (required)
    charge, rt, im                 (optional)

The A/B spike-in ratios are carried implicitly in the reported intensities and
organism membership — there is NO proteome_mix / dilution step. Run once per
sample (Sample A, Sample B).

Ute's DIA-NN reports are long-format: BOTH conditions (runs `A1_*` and `B1_*`)
and all 6 replicates are concatenated in a single file. Use `--run-prefix` to
select one condition's runs (defaults to `<sample>1`, i.e. A->A1 / B->B1);
intensities are then taken as the per-precursor MEDIAN across that condition's
replicates.

Key contract details (see plasmabench/peptides.py and TimSim load_findings):
  - `sequence` must be UniMod-notated and parseable by the Rust PeptideSequence
    parser; rows that don't validate are dropped (with a count) so the Rust core
    never aborts.
  - Optional columns (rt/im) are written as real values or the column is OMITTED
    entirely — never blank strings. TimSim runs dropna() over every present
    column, so one all-empty optional column would delete every row.

Usage (single combined DIA-NN report, selecting each condition by run prefix):
    python scripts/build_seed_from_report.py \
        --engine diann --report .../G1_nLC_tTOF_report.tsv \
        --sample A --run-prefix A1 \
        --out simulations/stage1/seeds/seed_sampleA.csv

    python scripts/build_seed_from_report.py \
        --engine diann --report .../G1_nLC_tTOF_report.tsv \
        --sample B --run-prefix B1 \
        --out simulations/stage1/seeds/seed_sampleB.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.peptides import (  # noqa: E402
    diann_modseq_to_unimod,
    fragpipe_to_unimod,
    is_valid_timsim_sequence,
)
from plasmabench.species import assign_species, PYE1_COMPOSITION  # noqa: E402
from plasmabench.intensities import baseline_intensity, get_tenzer_curve  # noqa: E402

# Column name candidates per engine. Search-engine versions differ; the first
# present candidate wins, and we error clearly if none match (so a schema drift
# fails loudly instead of silently dropping data). Validate against Ute's
# actual exports — see TODO.md.
DIANN_COLS = {
    # Species tags (`_HUMAN`/`_YEAS8`/`_ECOLI`) live in Protein.Names, NOT in the
    # bare accessions of Protein.Group — so the seed's `protein` field is built
    # from Names (falling back to the accession row-wise where Names is empty,
    # e.g. FragPipe's DIA-NN-format export leaves Names blank).
    "protein_names": ["Protein.Names"],
    "protein_acc": ["Protein.Group", "Protein.Ids"],
    "modseq": ["Modified.Sequence"],
    "charge": ["Precursor.Charge"],
    "intensity": ["Precursor.Quantity", "Precursor.Normalised", "Ms1.Area"],
    "rt_min": ["RT", "RT.Start"],          # DIA-NN reports RT in MINUTES
    "im": ["IM", "iIM", "Predicted.IM"],
    "run": ["Run"],                        # per-injection label, e.g. A1_G_DIA_..._R3
}
FRAGPIPE_COLS = {
    "protein": ["Protein", "Protein ID"],
    "peptide": ["Peptide", "Peptide Sequence"],
    "mods": ["Assigned Modifications"],
    "charge": ["Charge"],
    "intensity": ["Intensity"],
    "rt_sec": ["Retention"],               # FragPipe psm.tsv Retention is SECONDS
    "im": ["Ion Mobility"],
    "run": ["Run", "Spectrum File"],
}


def _pick(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _read(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".parquet", ".pq"):
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def load_diann(df: pd.DataFrame) -> pd.DataFrame:
    cols = {k: _pick(df, v) for k, v in DIANN_COLS.items()}
    for required in ("modseq", "intensity"):
        if cols[required] is None:
            raise SystemExit(f"DIA-NN report missing a column for '{required}' "
                             f"(looked for {DIANN_COLS[required]})")
    if cols["protein_names"] is None and cols["protein_acc"] is None:
        raise SystemExit("DIA-NN report missing a protein column "
                         f"(looked for {DIANN_COLS['protein_names']} / "
                         f"{DIANN_COLS['protein_acc']})")
    # protein = species-tagged Protein.Names, row-wise falling back to the
    # accession where Names is blank (keeps the species oracle fed; see
    # plasmabench/species.py).
    names = (df[cols["protein_names"]].astype(str).str.strip()
             if cols["protein_names"] else pd.Series("", index=df.index))
    acc = (df[cols["protein_acc"]].astype(str).str.strip()
           if cols["protein_acc"] else pd.Series("", index=df.index))
    blank = names.isin(("", "nan", "NaN", "None"))
    protein = names.mask(blank, acc)
    out = pd.DataFrame({
        "protein": protein,
        "sequence": df[cols["modseq"]].astype(str).map(diann_modseq_to_unimod),
        "intensity": pd.to_numeric(df[cols["intensity"]], errors="coerce"),
    })
    if cols["charge"]:
        out["charge"] = pd.to_numeric(df[cols["charge"]], errors="coerce")
    if cols["rt_min"]:
        out["rt"] = pd.to_numeric(df[cols["rt_min"]], errors="coerce") * 60.0  # min → s
    if cols["im"]:
        out["im"] = pd.to_numeric(df[cols["im"]], errors="coerce")
    return out


def load_fragpipe(df: pd.DataFrame) -> pd.DataFrame:
    cols = {k: _pick(df, v) for k, v in FRAGPIPE_COLS.items()}
    for required in ("protein", "peptide", "intensity"):
        if cols[required] is None:
            raise SystemExit(f"FragPipe report missing a column for '{required}' "
                             f"(looked for {FRAGPIPE_COLS[required]})")
    mods = df[cols["mods"]] if cols["mods"] else pd.Series([""] * len(df))
    out = pd.DataFrame({
        "protein": df[cols["protein"]].astype(str),
        "sequence": [fragpipe_to_unimod(p, m)
                     for p, m in zip(df[cols["peptide"]].astype(str), mods.fillna(""))],
        "intensity": pd.to_numeric(df[cols["intensity"]], errors="coerce"),
    })
    if cols["charge"]:
        out["charge"] = pd.to_numeric(df[cols["charge"]], errors="coerce")
    if cols["rt_sec"]:
        out["rt"] = pd.to_numeric(df[cols["rt_sec"]], errors="coerce")
    if cols["im"]:
        out["im"] = pd.to_numeric(df[cols["im"]], errors="coerce")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True, choices=["diann", "fragpipe"])
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sample", required=True, choices=["A", "B"],
                    help="PYE1 sample label — drives the species-audit expectation "
                         "and the default --run-prefix")
    ap.add_argument("--run-prefix", default=None,
                    help="Keep only rows whose Run label STARTS WITH this prefix "
                         "(e.g. 'A1' or 'B1'). Ute's reports concatenate BOTH "
                         "conditions x 6 replicates in one file; without this the "
                         "seed merges A and B and destroys the ratio. Defaults to "
                         "'<sample>1' (A->A1, B->B1); pass '' to disable.")
    ap.add_argument("--run-contains", default=None,
                    help="Keep only rows whose Run label CONTAINS this substring "
                         "(takes precedence over --run-prefix). Use for the L-site "
                         "reports whose runs are slot-coded mid-string, e.g. "
                         "'Slot2-11' (=Sample A) / 'Slot2-12' (=Sample B).")
    ap.add_argument("--keep-species", nargs="+", default=None,
                    choices=["HUMAN", "YEAST", "ECOLI"],
                    help="Keep only peptides assigned to these species; everything "
                         "else (other species, plus ambiguous/contaminant rows that "
                         "resolve to None) is dropped. Use '--keep-species YEAST ECOLI' "
                         "to build YE-only seeds for the Stage-2 superimpose-on-real-"
                         "plasma route, where the human background comes from the real "
                         "plasma .d rather than the seed.")
    ap.add_argument("--ratio-mode", choices=["observed", "imposed"], default="observed",
                    help="observed: use the report's (DIA-NN-compressed) intensities; "
                         "the oracle is the seeded ratio. imposed: discard found "
                         "intensities and assign baseline x species dilution factor "
                         "(the HYE route) so the A/B ratio is exactly nominal "
                         "(human 1.0 / yeast 3.0 / ecoli 0.5) — a clean accuracy oracle.")
    ap.add_argument("--rt-max", type=float, default=None,
                    help="Linearly rescale the rt column so its max maps to this "
                         "value (seconds). Use to fit seeds onto a SHORTER reference "
                         "gradient: the G-site IDs elute to ~2121s (35min) but the "
                         "L-method DIA blank spans only ~1500s (25min), so without "
                         "this TimSim drops ~26%% of peptides (rt outside the "
                         "reference frame range). RT is not a scored quantity, so a "
                         "uniform compression is safe and keeps all IDs.")
    ap.add_argument("--upscale-factor", type=float, default=1e5,
                    help="imposed mode: Tenzer baseline scale (matches TimSim HYE "
                         "default 1e5). Absolute scale is re-normalized by TimSim's "
                         "from_findings loader, so this is mostly cosmetic.")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    # --run-contains (substring) wins over --run-prefix (startswith). Only fall
    # back to the '<sample>1' prefix default when neither is given.
    run_prefix = (f"{args.sample}1" if args.run_prefix is None else args.run_prefix) \
        if args.run_contains is None else None

    df = _read(args.report)

    # Filter to one condition's runs. These reports are long-format with both
    # conditions x 6 replicates concatenated (A1_*/B1_* at site G, Slot2-11/
    # Slot2-12 at site L).
    run_col = _pick(df, DIANN_COLS["run"] if args.engine == "diann"
                    else FRAGPIPE_COLS["run"])
    if args.run_contains:
        if run_col is None:
            raise SystemExit(f"--run-contains '{args.run_contains}' given but no Run "
                             f"column found (engine={args.engine}).")
        runs = df[run_col].astype(str)
        keep = runs.str.contains(args.run_contains, regex=False)
        if not keep.any():
            raise SystemExit(f"--run-contains '{args.run_contains}' matched 0 of "
                             f"{len(df):,} rows. Distinct runs: {sorted(runs.unique())[:20]}")
        kept_runs = sorted(runs[keep].unique())
        df = df[keep]
        print(f"  run filter '*{args.run_contains}*' kept {len(df):,} rows "
              f"across {len(kept_runs)} run(s): {kept_runs}")
    elif run_prefix:
        if run_col is None:
            raise SystemExit(f"--run-prefix '{run_prefix}' given but no Run column "
                             f"found to filter on (engine={args.engine}).")
        runs = df[run_col].astype(str)
        keep = runs.str.startswith(run_prefix)
        if not keep.any():
            raise SystemExit(f"--run-prefix '{run_prefix}' matched 0 of {len(df):,} "
                             f"rows. Distinct runs: {sorted(runs.unique())[:20]}")
        kept_runs = sorted(runs[keep].unique())
        df = df[keep]
        print(f"  run filter '{run_prefix}*' kept {len(df):,} rows "
              f"across {len(kept_runs)} run(s): {kept_runs}")
    elif run_col is not None and df[run_col].astype(str).nunique() > 1:
        print(f"  WARNING: no --run-prefix set but report has "
              f"{df[run_col].nunique()} distinct runs — A and B may be merged.")

    out = load_diann(df) if args.engine == "diann" else load_fragpipe(df)
    n_raw = len(out)

    # Always drop rows with a sequence TimSim can't parse. The intensity>0 filter
    # is a detectability criterion that only applies to OBSERVED mode (which uses
    # the reported intensity); imposed mode discards it, so an identified
    # precursor with a missing/zero reported quantity must still be kept.
    if args.ratio_mode == "observed":
        out = out[out["intensity"].notna() & (out["intensity"] > 0)]
    bad_seq = ~out["sequence"].map(is_valid_timsim_sequence)
    n_bad = int(bad_seq.sum())
    out = out[~bad_seq]

    # Aggregate to one row per (sequence, charge). After the run filter each key
    # appears at most once per replicate; take the MEDIAN intensity across
    # replicates (robust, and unlike a sum it isn't inflated by how many
    # replicates happened to detect the precursor — which would bias the A/B
    # ratio toward more-frequently-observed peptides). rt/im likewise median.
    # `protein` is the UNION of every protein string seen for the key (not the
    # first row): a replicate-order-dependent pick could otherwise assign the
    # wrong species, and differently in A vs B. assign_species() then resolves
    # cross-species unions to None (dropped from the oracle).
    def _union(s):
        return ";".join(dict.fromkeys(str(x) for x in s if pd.notna(x) and str(x)))
    has_charge = "charge" in out.columns and out["charge"].notna().any()
    group_keys = ["sequence", "charge"] if has_charge else ["sequence"]
    agg = {"intensity": ("intensity", "median"),
           "protein": ("protein", _union)}
    if "rt" in out.columns:
        agg["rt"] = ("rt", "median")
    if "im" in out.columns:
        agg["im"] = ("im", "median")
    # dropna=False so a stray NaN charge becomes its own group rather than being
    # silently deleted (the optional-column logic below then omits charge).
    grouped = out.groupby(group_keys, as_index=False, dropna=False).agg(**agg)

    # Optional species restriction (applies to BOTH ratio modes). The union
    # `protein` is resolved by assign_species(), so cross-species peptides (→ None)
    # are dropped here too. YE-only seeds for Stage 2 keep {YEAST, ECOLI}.
    n_species_dropped = 0
    if args.keep_species:
        keep_set = set(args.keep_species)
        keep_mask = grouped["protein"].map(assign_species).isin(keep_set)
        n_species_dropped = int((~keep_mask).sum())
        grouped = grouped[keep_mask].reset_index(drop=True)
        if grouped.empty:
            raise SystemExit(f"--keep-species {sorted(keep_set)} matched 0 rows "
                             f"(dropped all {n_species_dropped:,}). Wrong species set?")

    # --- imposed ratio mode -------------------------------------------------
    # Keep the found IDs (sequence/charge/rt/im) but THROW AWAY the found
    # intensity. Each ion gets a Tenzer baseline (deterministic per sequence+
    # charge, so A and B share it) scaled by the per-species mass fraction for
    # this sample. The A/B fold change is then exactly PYE1's nominal ratio,
    # free of the search engine's quant compression. Mirrors TimSim's HYE
    # proteome_mix route; see plasmabench/intensities.py.
    n_imposed_dropped = 0
    if args.ratio_mode == "imposed":
        sp = grouped["protein"].map(assign_species)
        keep = sp.notna()                       # drop ambiguous / contaminants
        n_imposed_dropped = int((~keep).sum())
        grouped = grouped[keep].reset_index(drop=True)
        sp = sp[keep].reset_index(drop=True)
        comp = PYE1_COMPOSITION[args.sample]    # HUMAN/YEAST/ECOLI mass fractions
        curve = get_tenzer_curve()
        charges = (grouped["charge"] if "charge" in grouped.columns
                   else pd.Series([None] * len(grouped)))
        grouped["intensity"] = [
            baseline_intensity(s, c, upscale_factor=args.upscale_factor, curve=curve)
            * comp[org]
            for s, c, org in zip(grouped["sequence"], charges, sp)
        ]

    # Optionally rescale RT onto a shorter reference gradient (see --rt-max).
    rt_rescale = None
    if args.rt_max is not None and "rt" in grouped.columns and len(grouped):
        rt_observed_max = float(grouped["rt"].max())
        if rt_observed_max > 0:
            rt_rescale = args.rt_max / rt_observed_max
            grouped["rt"] = grouped["rt"] * rt_rescale

    # Assemble in contract order; only include optional columns that are fully
    # populated — a partially-NaN optional column would drop rows in TimSim.
    seed = pd.DataFrame({
        "protein": grouped["protein"].fillna(""),
        "sequence": grouped["sequence"],
        "intensity": grouped["intensity"].astype(float),
    })
    dropped_optional = []
    for opt, caster in (("charge", lambda s: s.astype("Int64")),
                        ("rt", lambda s: s.astype(float)),
                        ("im", lambda s: s.astype(float))):
        if opt in grouped.columns and grouped[opt].notna().all():
            seed[opt] = caster(grouped[opt])
        elif opt in grouped.columns:
            dropped_optional.append(opt)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    seed.to_csv(args.out, index=False)

    # Species audit — the ratio oracle depends on a clean organism split.
    species = grouped["protein"].map(assign_species)
    counts = species.value_counts(dropna=False).to_dict()
    print(f"wrote {args.out}")
    print(f"  engine={args.engine} sample={args.sample} ratio-mode={args.ratio_mode}")
    print(f"  {len(seed):,} rows  ({grouped['sequence'].nunique():,} unique sequences)"
          f"  from {n_raw:,} raw report rows")
    if n_bad:
        print(f"  dropped {n_bad:,} rows with sequences TimSim can't parse")
    if args.keep_species:
        print(f"  --keep-species {sorted(set(args.keep_species))}: dropped "
              f"{n_species_dropped:,} rows of other/ambiguous species")
    if args.ratio_mode == "imposed":
        print(f"  imposed mode: intensity = Tenzer baseline x species fraction "
              f"(upscale={args.upscale_factor:g}); found intensities discarded")
        if n_imposed_dropped:
            print(f"  dropped {n_imposed_dropped:,} rows with no clean species "
                  f"(ambiguous / contaminant) — excluded from imposed truth")
    if rt_rescale is not None:
        print(f"  rescaled rt by x{rt_rescale:.4f} -> max {grouped['rt'].max():.0f}s "
              f"(fit onto shorter reference gradient)")
    if dropped_optional:
        print(f"  omitted optional column(s) {dropped_optional} (not fully populated)")
    print(f"  species split: {counts}")
    print("  expected PYE1 mass fractions for sample "
          f"{args.sample}: HUMAN .90 / "
          f"{'YEAST .02 / ECOLI .08' if args.sample == 'A' else 'YEAST .06 / ECOLI .04'}")


if __name__ == "__main__":
    main()
