Reading additional input from stdin...
OpenAI Codex v0.136.0
--------
workdir: /scratch/timsim-demo/PlasmaBENCH
model: gpt-5.5
provider: openai
approval: never
sandbox: read-only
reasoning effort: medium
reasoning summaries: none
session id: 019ea7c7-ea69-7c42-8f35-e918945c5971
--------
user
You are an independent reviewer of PlasmaBENCH, a benchmark that seeds the TimSim mass-spec simulator from real DIA-NN search reports (the 'from_findings' mode) to build ground-truth DIA-PASEF runs for a 3-species (human/yeast/ecoli) plasma quant benchmark (PYE1: A vs B with nominal species ratios human 1.0, yeast 3.0, ecoli 0.5).

Review the bundled source files for CORRECTNESS and robustness. Key context to judge against:
- build_seed_from_report.py emits a 6-col findings CSV (protein,sequence,intensity[,charge,rt,im]). Two ratio modes: 'observed' keeps the report's (DIA-NN-compressed) intensities; 'imposed' DISCARDS found intensity and assigns intensity = Tenzer-baseline (deterministic per sequence+charge via crc32 so Sample A and B share it) x per-species mass fraction, to produce an exactly-nominal A/B ratio. This mirrors TimSim's HYE proteome_mix route.
- IMPORTANT subtlety I want validated: TimSim's load_findings normalizes each seed file by its median sequence-intensity, then x upscale_factor. Since A and B are separate files, a global offset can appear; we claim human-reference (LFQbench-style) normalization removes it and recovers exact 1.0/3.0/0.5. Is that reasoning sound? Any case where imposed mode does NOT give exact ratios?
- species.py assigns organism by UniProt entry-name suffix; yeast is _YEAS8 (EC1118 strain), not _YEAST.
- get_fasta.sh reconstructs the search DB from UniProt; qc_fasta_coverage.py digests it (Trypsin/P) to measure seed-peptide coverage.

Focus on: 1) seed-building correctness (run filtering, charge aggregation, dropna/optional-column handling, the imposed baseline+fraction math, determinism across A/B), 2) species/intensity math edge cases, 3) the digest in qc_fasta_coverage.py (is it a faithful Trypsin/P digest? off-by-one in cleavage sites? does it match how DIA-NN would search?), 4) reproducibility-script robustness (get_fasta.sh / get_diann.sh / run.sh: error handling, silent failures). Give concrete, specific, actionable findings with file:concept references. Flag real bugs prominently. Cap at ~900 words.

<stdin>


===================== FILE: scripts/build_seed_from_report.py =====================
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

    # Drop rows without the required fields (intensity > 0, real sequence).
    out = out[out["intensity"].notna() & (out["intensity"] > 0)]
    bad_seq = ~out["sequence"].map(is_valid_timsim_sequence)
    n_bad = int(bad_seq.sum())
    out = out[~bad_seq]

    # Aggregate to one row per (sequence, charge). After the run filter each key
    # appears at most once per replicate; take the MEDIAN intensity across
    # replicates (robust, and unlike a sum it isn't inflated by how many
    # replicates happened to detect the precursor — which would bias the A/B
    # ratio toward more-frequently-observed peptides). rt/im likewise median.
    has_charge = "charge" in out.columns and out["charge"].notna().any()
    group_keys = ["sequence", "charge"] if has_charge else ["sequence"]
    agg = {"intensity": ("intensity", "median"),
           "protein": ("protein", lambda s: s.iloc[0])}
    if "rt" in out.columns:
        agg["rt"] = ("rt", "median")
    if "im" in out.columns:
        agg["im"] = ("im", "median")
    grouped = out.groupby(group_keys, as_index=False).agg(**agg)

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


===================== FILE: plasmabench/species.py =====================
"""Organism assignment and the PYE1 expected-ratio table.

The three-species design lives or dies on a clean organism map: every
quantified peptide must be attributable to human (plasma), yeast, or
E. coli, and the few that aren't (shared/ambiguous, contaminants) must be
excluded from the ratio oracle rather than silently miscounted.

Assignment is by UniProt entry-name suffix (`_HUMAN`, `_YEAS8`, `_ECOLI`),
which is what DIA-NN and FragPipe carry through in their `Protein.Names`
column when the search FASTA is the standard SwissProt download per
organism. If Ute's FASTA uses different organism tags, extend `SPECIES_TAGS`.

NOTE: tags live in `Protein.Names` (e.g. `ALBU_HUMAN`), NOT `Protein.Group`
(bare accessions like `P02768`). Feed the names column here.

The yeast proteome Ute used is the *strain-specific* S. cerevisiae download,
whose UniProt OS code is `YEAS8` (strain ATCC 204508 / S288c) — so the tag is
`_YEAS8`, not the canonical `_YEAST`. Both are mapped to YEAST so a mixed
FASTA still resolves. Contaminant entries (`_CONTA`, the cRAP set) match no
tag and therefore return None — correctly excluded from the ratio oracle.
"""
from __future__ import annotations

import re

# UniProt OS tags → canonical species label. E. coli SwissProt entries are
# `*_ECOLI`; human `*_HUMAN`; yeast appears as `*_YEAS8` (strain-specific) or
# `*_YEAST` (canonical) — both map to YEAST.
SPECIES_TAGS: dict[str, str] = {
    "_HUMAN": "HUMAN",
    "_YEAS8": "YEAST",
    "_YEAST": "YEAST",
    "_ECOLI": "ECOLI",
}

# PYE1 nominal composition (fraction of total protein mass per sample).
PYE1_COMPOSITION: dict[str, dict[str, float]] = {
    "A": {"HUMAN": 0.90, "YEAST": 0.02, "ECOLI": 0.08},
    "B": {"HUMAN": 0.90, "YEAST": 0.06, "ECOLI": 0.04},
}

# Expected Sample-B / Sample-A intensity ratio per species — the quant oracle.
EXPECTED_RATIO_B_OVER_A: dict[str, float] = {
    sp: PYE1_COMPOSITION["B"][sp] / PYE1_COMPOSITION["A"][sp]
    for sp in PYE1_COMPOSITION["A"]
}  # HUMAN 1.0, YEAST 3.0, ECOLI 0.5


def assign_species(protein_field: str) -> str | None:
    """Return 'HUMAN' | 'YEAST' | 'ECOLI', or None if absent/ambiguous.

    A peptide mapping to more than one organism (protein groups that span
    species) is ambiguous and returns None — such peptides must not enter
    the ratio oracle.

    Examples:
        >>> assign_species("sp|P02768|ALBU_HUMAN")
        'HUMAN'
        >>> assign_species("P00330_YEAST;P00331_YEAST")
        'YEAST'
        >>> assign_species("sp|P02768|ALBU_HUMAN;sp|P0A6F5|CH60_ECOLI") is None
        True
    """
    if not isinstance(protein_field, str) or not protein_field:
        return None
    found: set[str] = set()
    for tag, label in SPECIES_TAGS.items():
        if tag in protein_field:
            found.add(label)
    if len(found) == 1:
        return next(iter(found))
    return None  # zero matches or cross-species → ambiguous


def expected_log2_ratio(species: str) -> float | None:
    """log2(Sample B / Sample A) expected for a species, or None if unknown."""
    import math
    r = EXPECTED_RATIO_B_OVER_A.get(species)
    return None if r is None else math.log2(r)


__all__ = [
    "SPECIES_TAGS",
    "PYE1_COMPOSITION",
    "EXPECTED_RATIO_B_OVER_A",
    "assign_species",
    "expected_log2_ratio",
]


===================== FILE: plasmabench/intensities.py =====================
"""Synthetic per-peptide baseline intensities for the `imposed` ratio mode.

The `imposed` seed reproduces what TimSim's HYE / proteome_mix pipeline does for
intensity: it does NOT use the search engine's reported (and ratio-compressed)
intensity. Instead each peptide draws a baseline abundance from the Tenzer
rank-abundance curve, and the A/B fold change comes purely from a per-species
dilution factor (see plasmabench/species.py PYE1_COMPOSITION). The result is a
clean, uncompressed ground truth.

Ported verbatim from rustims imspy_simulation/timsim/jobs/simulate_proteins.py
(`get_tenzer_hokey` + `assign_events`) so the two stay in sync:

    events = upscale_factor * tenzer_curve[random_rank]

The one change for seeding is determinism: the random rank is drawn from a
per-peptide RNG keyed on (sequence, charge) via a stable hash, so an independent
Sample-A and Sample-B invocation assign the SAME baseline to a shared peptide —
which is what makes the imposed fold change exact (mirrors HYE's
`is_fold_change_condition`, where condition B reuses condition A's baseline).
"""
from __future__ import annotations

import zlib

import numpy as np

_CURVE_LEN = int(1e4)


def get_tenzer_curve() -> np.ndarray:
    """The Tenzer rank-abundance curve (length 1e4), highest abundance at rank 1.

    Verbatim from rustims `get_tenzer_hokey()`.
    """
    delta_exp = 0.005
    delta = 0.0011
    exp_red = 600
    target_start = 1e4

    rank = np.linspace(1, _CURVE_LEN, _CURVE_LEN)
    exp_values = 1 / np.exp(rank / exp_red)
    cumulative_exp_sum = np.cumsum(exp_values)
    return target_start * (2 ** (-delta * rank)) * (2 ** (-delta_exp * cumulative_exp_sum))


def _stable_seed(sequence: str, charge) -> int:
    """Deterministic 32-bit seed for a peptide ion, stable across processes.

    Python's built-in hash() is salted per-process, so it would give Sample A
    and Sample B different baselines — use crc32 instead.
    """
    return zlib.crc32(f"{sequence}|{charge}".encode("utf-8")) & 0xFFFFFFFF


def baseline_intensity(sequence: str, charge, *, upscale_factor: float,
                       curve: np.ndarray | None = None) -> float:
    """Per-peptide baseline abundance = upscale_factor * tenzer[random_rank].

    The rank is drawn deterministically from (sequence, charge) so the same ion
    gets the same baseline in every sample — the basis of an exact fold change.
    """
    if curve is None:
        curve = get_tenzer_curve()
    rng = np.random.RandomState(_stable_seed(sequence, charge))
    rank = int(rng.uniform(1, _CURVE_LEN))
    rank = min(rank, _CURVE_LEN - 1)
    return float(upscale_factor) * float(curve[rank])


__all__ = ["get_tenzer_curve", "baseline_intensity"]


===================== FILE: scripts/render_configs.py =====================
"""Render per-sample TimSim configs from a base template.

The base config uses container-shaped paths:
    /work = repo root
    /data = ${PLB_DATA:-<repo>/data}        (raw blanks + plasma + fastas)
    /sim  = ${PLB_OUT:-<repo>/simulations}  (TimSim outputs)
That works inside the rustims Docker image (volumes mounted by Makefile /
docker-compose). From a host venv (`make venv-source`) those paths don't
exist, so `--path-mode host` rewrites them to host equivalents.

Substitutes SAMPLE_PLACEHOLDER and repoints findings_path / save_path /
experiment_name for each sample.

Usage:
    python scripts/render_configs.py \
        --base simulations/stage1/configs/base-dia.toml \
        --seeds-dir simulations/stage1/seeds \
        --out-dir simulations/stage1/dia \
        --samples A B \
        --path-mode host
"""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


def host_substitutions(repo_root: Path) -> dict[str, str]:
    repo = repo_root.resolve()
    data = os.environ.get("PLB_DATA") or str(repo / "data")
    sim = os.environ.get("PLB_OUT") or str(repo / "simulations")
    return {"/work": str(repo), "/data": data, "/sim": sim}


def rewrite_paths(text: str, subs: dict[str, str]) -> str:
    out = text
    for src, dst in subs.items():
        pattern = re.compile(rf'"({re.escape(src)})((?:/[^"]*)?)"')
        out = pattern.sub(lambda m: f'"{dst}{m.group(2)}"', out)
    return out


def render(base_toml: str, *, sample: str, seed_csv: str, save_path: str,
           experiment_name: str) -> str:
    out = base_toml.replace("SAMPLE_PLACEHOLDER", sample)
    new_lines = []
    for line in out.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("findings_path"):
            new_lines.append(f'findings_path   = "{seed_csv}"')
        elif stripped.startswith("save_path"):
            new_lines.append(f'save_path       = "{save_path}"')
        elif stripped.startswith("experiment_name"):
            new_lines.append(f'experiment_name      = "{experiment_name}"')
        else:
            new_lines.append(line)
    return "\n".join(new_lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--seeds-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--samples", nargs="+", default=["A", "B"])
    ap.add_argument("--seed-suffix", default="",
                    help="Seed/variant suffix, e.g. '_imposed' to use "
                         "seed_sample<A|B>_imposed.csv and write to sample<A|B>_imposed/ "
                         "(keeps observed and imposed runs from colliding).")
    ap.add_argument("--path-mode", choices=["container", "host"], default="container")
    ap.add_argument("--repo-root", type=Path, default=None)
    ap.add_argument("--container-base", default="/work")
    args = ap.parse_args()

    repo_root = args.repo_root or Path(__file__).resolve().parents[1]
    base_text = args.base.read_text(encoding="utf-8")

    if args.path_mode == "host":
        subs = host_substitutions(repo_root)
        base_text = rewrite_paths(base_text, subs)
        work_prefix = subs["/work"]
    else:
        work_prefix = args.container_base

    rel_seed = args.seeds_dir.relative_to(repo_root) \
        if args.seeds_dir.is_absolute() else args.seeds_dir
    rel_save = args.out_dir.relative_to(repo_root) \
        if args.out_dir.is_absolute() else args.out_dir

    suffix = args.seed_suffix
    for sample in args.samples:
        sample_dir = args.out_dir / f"sample{sample}{suffix}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        seed_path = f"{work_prefix}/{rel_seed}/seed_sample{sample}{suffix}.csv"
        save_path = f"{work_prefix}/{rel_save}/sample{sample}{suffix}"
        experiment = f"PLB-S1-DIA-sample{sample}{suffix.replace('_', '-')}"
        text = render(base_text, sample=sample, seed_csv=seed_path,
                      save_path=save_path, experiment_name=experiment)
        out_path = sample_dir / "config.toml"
        out_path.write_text(text, encoding="utf-8")
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()


===================== FILE: configs/dia-nn/diann-pye1.cfg =====================
--fasta-search
--gen-spec-lib
--predictor
--cut K*,R*
--missed-cleavages 1
--min-pep-len 7
--max-pep-len 30
--min-pr-charge 1
--max-pr-charge 4
--min-pr-mz 300
--max-pr-mz 1800
--min-fr-mz 200
--max-fr-mz 1800
--unimod4
--var-mods 0
--qvalue 0.01
--matrix-qvalue 0.01
--reanalyse
--mass-acc 0
--mass-acc-ms1 0
--window 0
--smart-profiling


===================== FILE: configs/dia-nn/run.sh =====================
#!/usr/bin/env bash
# Run a DIA-NN version on PlasmaBENCH SIM .d files with the UNIFIED PYE1 search
# settings (configs/dia-nn/diann-pye1.cfg). Both versions get the identical
# settings file — the only intended differences are the binary and its built-in
# defaults (QuantUMS, DL models): that IS the version comparison.
#
# Pass Sample A and Sample B .d together so DIA-NN's cross-run / MBR step
# (--reanalyse, in the cfg) can build the A-vs-B quant matrix.
#
# Usage:
#   bash configs/dia-nn/run.sh <1.8|2.5> <fasta> <out_dir> <runA.d> <runB.d> [...]
# Env:
#   THREADS   (default: nproc)
#
# Example (once the SIM .d's and FASTA exist):
#   bash configs/dia-nn/run.sh 1.8 data/raw/fasta/plasma_yeast_ecoli.fasta \
#     results/diann-1.8 \
#     simulations/stage1/dia/sampleA_imposed/*.d \
#     simulations/stage1/dia/sampleB_imposed/*.d
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)

[[ $# -ge 4 ]] || { echo "usage: run.sh <1.8|2.5> <fasta> <out_dir> <d1> [d2 ...]" >&2; exit 1; }
VER="$1"; FASTA="$2"; OUT_DIR="$3"; shift 3

case "$VER" in
  1.8|1.8.1) BIN="$REPO_ROOT/tools/diann/1.8.1/diann-1.8.1" ;;
  2.5|2.5.0) BIN="$REPO_ROOT/tools/diann/2.5.0/diann-2.5.0/diann-linux" ;;
  *) echo "unknown DIA-NN version '$VER' (use 1.8 or 2.5)" >&2; exit 1 ;;
esac
[[ -x "$BIN" ]] || { echo "missing/!exec $BIN — run: bash scripts/get_diann.sh" >&2; exit 2; }
[[ -e "$FASTA" ]] || { echo "FASTA not found: $FASTA (TODO: plasma+yeast+ecoli reference)" >&2; exit 2; }

CFG="$REPO_ROOT/configs/dia-nn/diann-pye1.cfg"
SETTINGS=$(grep -vE '^\s*(#|$)' "$CFG" | tr '\n' ' ')
THREADS="${THREADS:-$(nproc)}"
mkdir -p "$OUT_DIR"

F_ARGS=(); for d in "$@"; do F_ARGS+=(--f "$d"); done

echo "==> DIA-NN $VER  ($BIN)"
echo "    fasta=$FASTA  out=$OUT_DIR  threads=$THREADS  runs=$#"
set -x
"$BIN" $SETTINGS \
  --fasta "$FASTA" \
  "${F_ARGS[@]}" \
  --out "$OUT_DIR/report.tsv" \
  --out-lib "$OUT_DIR/report-lib" \
  --temp "$OUT_DIR" \
  --threads "$THREADS"


===================== FILE: scripts/get_diann.sh =====================
#!/usr/bin/env bash
# Fetch the DIA-NN Linux builds used as software-under-test (versions 1.8 and 2.5).
#
# Binaries are NOT committed (tools/ is gitignored); this script makes the
# download reproducible. Academic builds from the official GitHub releases
# (https://github.com/vdemichev/DiaNN). Review the bundled LICENSE.txt before use.
#
# Usage:  bash scripts/get_diann.sh
# Result: tools/diann/1.8.1/  and  tools/diann/2.5.0/  with the diann binary inside.
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
DEST="$REPO_ROOT/tools/diann"
mkdir -p "$DEST"

# version label -> release asset URL
declare -A URLS=(
  ["1.8.1"]="https://github.com/vdemichev/DiaNN/releases/download/1.8.1/diann_1.8.1.tar.gz"
  ["2.5.0"]="https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.5.0-Academia-Linux.zip"
)

for ver in "${!URLS[@]}"; do
  url="${URLS[$ver]}"
  out="$DEST/$ver"
  if [[ -e "$out/.ok" ]]; then echo "==> $ver already present, skipping"; continue; fi
  mkdir -p "$out"
  archive="$DEST/$(basename "$url")"
  echo "==> downloading DIA-NN $ver"
  curl -L --fail --retry 3 -o "$archive" "$url"
  echo "==> extracting into $out"
  case "$archive" in
    *.tar.gz) tar -xzf "$archive" -C "$out" ;;
    *.zip)    unzip -q -o "$archive" -d "$out" ;;
  esac
  rm -f "$archive"
  touch "$out/.ok"
done

echo
echo "Installed:"
find "$DEST" -maxdepth 3 -type f \( -name 'diann*' -o -name '*.exe' \) -perm -u+x 2>/dev/null | sort || true
find "$DEST" -maxdepth 3 -iname 'diann*' | sort


===================== FILE: scripts/get_fasta.sh =====================
#!/usr/bin/env bash
# Build the PlasmaBENCH search FASTA: human (plasma background) + S. cerevisiae
# EC1118 + E. coli K12 — the three proteomes behind Ute's PYE1 reports, picked to
# match the organism tags the seeds carry (_HUMAN / _YEAS8 / _ECOLI).
#
# Verified from the report accessions (scripts identify them via UniProt):
#   human   taxId 9606    reviewed Swiss-Prot   (P02768=ALBU_HUMAN, ...)
#   ecoli   taxId 83333   reviewed Swiss-Prot   K12, proteome UP000000625
#   yeast   taxId 643680  TrEMBL                 EC1118, proteome UP000000286 (_YEAS8)
#
# Output (gitignored): data/raw/fasta/plasma_yeast_ecoli.fasta
#
# NOTE: this is a best-effort reconstruction of Ute's search database. Prefer her
# EXACT FASTA if she can share it; run scripts/qc_fasta_coverage.py afterwards to
# measure how many seed peptides this database actually covers.
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR="$REPO_ROOT/data/raw/fasta"
mkdir -p "$OUT_DIR"
STREAM="https://rest.uniprot.org/uniprotkb/stream"

fetch() {  # <label> <query> <outfile>
  local label="$1" query="$2" out="$3"
  echo "==> $label"
  curl -sS --fail --retry 3 -G "$STREAM" \
    --data-urlencode "query=$query" \
    --data-urlencode "format=fasta" \
    --data-urlencode "includeIsoform=false" \
    -o "$out"
  echo "    $(grep -c '^>' "$out") sequences -> $out"
}

fetch "human (reviewed Swiss-Prot, 9606)" \
      "(organism_id:9606) AND (reviewed:true)" "$OUT_DIR/_human.fasta"
fetch "E. coli K12 (reviewed, UP000000625)" \
      "(proteome:UP000000625) AND (reviewed:true)" "$OUT_DIR/_ecoli.fasta"
# Use the FULL set of EC1118 (_YEAS8) entries by organism, NOT the redundancy-
# reduced reference proteome UP000000286 (4514): the full set is ~5984 and the
# extra ~1470 entries are needed to cover the seed yeast peptides (the reference
# proteome alone only covers ~77%). See scripts/qc_fasta_coverage.py.
fetch "S. cerevisiae EC1118 (organism 643680, all _YEAS8)" \
      "(organism_id:643680)" "$OUT_DIR/_yeast.fasta"

OUT="$OUT_DIR/plasma_yeast_ecoli.fasta"
cat "$OUT_DIR/_human.fasta" "$OUT_DIR/_yeast.fasta" "$OUT_DIR/_ecoli.fasta" > "$OUT"
rm -f "$OUT_DIR/_human.fasta" "$OUT_DIR/_ecoli.fasta" "$OUT_DIR/_yeast.fasta"

echo
echo "Combined: $OUT"
echo "  total sequences: $(grep -c '^>' "$OUT")"
echo "  by tag: HUMAN=$(grep -c '_HUMAN' "$OUT")  YEAS8=$(grep -c '_YEAS8' "$OUT")  ECOLI=$(grep -c '_ECOLI' "$OUT")"


===================== FILE: scripts/qc_fasta_coverage.py =====================
"""QC: does the search FASTA contain the seed peptides?

The seeds are the peptides TimSim will put into the SIM .d. If the FASTA
(searched by DIA-NN/FragPipe) doesn't contain a seed peptide's protein, that
peptide can never be identified — an ARTIFICIAL false negative that has nothing
to do with the software under test. This is the main risk of reconstructing
Ute's database from UniProt rather than using her exact FASTA (strain mismatch,
isoforms, contaminants). Run this to quantify the miss rate per species.

Method: in-silico Trypsin/P digest of the FASTA (missed cleavages 0..N, length
7..30), per species; then check what fraction of each seed's (mod-stripped)
sequences land in the matching species' peptide set.

Usage:
    python scripts/qc_fasta_coverage.py \
        --fasta data/raw/fasta/plasma_yeast_ecoli.fasta \
        --seeds simulations/stage1/seeds/seed_sampleA_imposed.csv \
                simulations/stage1/seeds/seed_sampleB_imposed.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.species import assign_species  # noqa: E402

_UNIMOD = re.compile(r"\[UNIMOD:\d+\]", re.IGNORECASE)


def strip_mods(seq: str) -> str:
    return _UNIMOD.sub("", str(seq))


def digest(seq: str, missed: int, min_len: int, max_len: int):
    """Trypsin/P: cut after every K and R. Yield peptides with 0..missed MC."""
    # cleavage sites = index after each K/R
    cuts = [0] + [i + 1 for i, a in enumerate(seq) if a in "KR"]
    if cuts[-1] != len(seq):
        cuts.append(len(seq))
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(cuts))):
            pep = seq[cuts[i]:cuts[j]]
            if min_len <= len(pep) <= max_len:
                yield pep


def build_peptide_sets(fasta: Path, missed: int, min_len: int, max_len: int):
    sets: dict[str, set] = {"HUMAN": set(), "YEAST": set(), "ECOLI": set()}
    header, sp, buf = None, None, []
    n_prot = 0

    def flush():
        nonlocal n_prot
        if sp in sets and buf:
            n_prot_local = "".join(buf)
            for pep in digest(n_prot_local, missed, min_len, max_len):
                sets[sp].add(pep)

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                sp = assign_species(header)
                buf = []
                n_prot += 1
            else:
                buf.append(line.strip())
    flush()
    return sets, n_prot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--seeds", required=True, nargs="+", type=Path)
    ap.add_argument("--missed-cleavages", type=int, default=2)
    ap.add_argument("--min-len", type=int, default=7)
    ap.add_argument("--max-len", type=int, default=30)
    args = ap.parse_args()

    print(f"digesting {args.fasta} (Trypsin/P, MC<= {args.missed_cleavages}, "
          f"len {args.min_len}-{args.max_len}) ...")
    sets, n_prot = build_peptide_sets(args.fasta, args.missed_cleavages,
                                      args.min_len, args.max_len)
    print(f"  {n_prot:,} proteins -> peptides per species: "
          + ", ".join(f"{k}={len(v):,}" for k, v in sets.items()))

    seed = pd.concat([pd.read_csv(s) for s in args.seeds], ignore_index=True)
    seed = seed.drop_duplicates(subset=["sequence"])
    seed["bare"] = seed["sequence"].map(strip_mods)
    seed["sp"] = seed["protein"].map(assign_species)

    print(f"\nseed peptides (unique sequences): {len(seed):,}")
    print(f"{'species':7} {'seed pep':>9} {'covered':>9} {'coverage':>9} {'missed':>8}")
    total_cov = total = 0
    for sp in ["HUMAN", "YEAST", "ECOLI"]:
        s = seed[seed["sp"] == sp]
        if not len(s):
            continue
        covered = s["bare"].isin(sets[sp]).sum()
        total_cov += covered
        total += len(s)
        miss = len(s) - covered
        print(f"{sp:7} {len(s):>9,} {covered:>9,} {100*covered/len(s):>8.1f}% {miss:>8,}")
    none = seed["sp"].isna().sum()
    print(f"\noverall coverage: {total_cov:,}/{total:,} ({100*total_cov/total:.1f}%)"
          f"   [{none:,} seed rows had no species tag, ignored]")
    print("\nMisses = seed peptides absent from the reconstructed DB -> would be "
          "artificial false negatives. High yeast miss rate => get Ute's exact "
          "EC1118 FASTA (strain variants).")


if __name__ == "__main__":
    main()
</stdin>
codex
**Findings**

1. **High: species assignment can become arbitrary during aggregation.**  
   `build_seed_from_report.py: groupby aggregation` selects `protein = s.iloc[0]`. If replicate rows report different protein groups, including cross-species assignments, the chosen species depends on row order. Imposed mode may therefore retain an ambiguous peptide and assign the wrong fraction, potentially differently in A and B. Aggregate all protein names for each `(sequence, charge)`, then call `assign_species` on the union; drop groups containing multiple species.

2. **High: partially missing charge silently deletes rows.**  
   `build_seed_from_report.py: has_charge/groupby` uses charge as a key whenever any charge is present. Pandas `groupby` drops rows whose charge is NaN, before the later optional-column logic can omit charge. Either require charge to be complete, use `dropna=False`, or explicitly separate/repair missing-charge rows. Also validate that charges are finite positive integers before casting to `Int64`.

3. **High: imposed mode is exact only on matched precursor ions, before independent seed normalization.**  
   For an ion present in both files with identical sequence, charge, and species assignment, the raw imposed ratio is exactly the composition ratio. TimSim’s independent median normalization changes it to:

   `observed B/A = nominal B/A × median(A seed intensities) / median(B seed intensities)`.

   A human-reference correction removes that common multiplier exactly in ideal arithmetic because every matched human ion has nominal ratio 1.0. The reasoning is therefore sound, but the claim should be qualified. Exact recovery fails or becomes empirical when:
   - A/B use different charge states or species assignments.
   - Normalization is estimated after missing values, stochastic simulation, censoring, saturation, or identification filtering.
   - There are too few shared human precursors.
   - “Human normalization” operates on protein summaries whose peptide sets differ.
   - Downstream processing applies nonlinear or run-specific normalization.

   Different A/B seed ID sets do not invalidate exactness for their intersection, but they can change both TimSim medians and the reference factor’s practical estimate.

4. **Medium: imposed mode unnecessarily depends on valid observed intensity.**  
   `build_seed_from_report.py: required-field filtering` drops rows with absent or nonpositive report intensity before imposed intensities are generated. This contradicts “discard found intensity”: an identified precursor with zero/missing quantity cannot enter imposed mode. In imposed mode, require identification fields but do not filter on reported intensity unless this is an intentional detectability criterion and is documented.

5. **Medium: sequence-level uniqueness hides charge-specific coverage and truth differences.**  
   `qc_fasta_coverage.py` and its reporting deduplicate only by `sequence`, while seeds and deterministic baselines are precursor-ion `(sequence, charge)` based. Peptide FASTA coverage is legitimately sequence-based, but report both unique peptides and affected seed ions; otherwise the apparent false-negative impact can be understated.

6. **Medium: QC digest default does not match the bundled DIA-NN search.**  
   `diann-pye1.cfg` specifies `--missed-cleavages 1`, while `qc_fasta_coverage.py` defaults to 2. This can label two-missed-cleavage seed peptides as searchable when DIA-NN will not search them. Change the default to 1 or pass search settings explicitly.

7. **Digest cleavage indexing itself is correct for Trypsin/P.**  
   `qc_fasta_coverage.py:digest` cuts after every K/R, including before proline, and the `j` bound correctly produces zero through `missed` missed cleavages. There is no off-by-one error. However, it is only an approximate DIA-NN search-space model: it does not validate amino-acid alphabet, terminal handling, precursor m/z/charge feasibility, modifications, or DIA-NN-specific peptide filtering. The configured `--cut K*,R*` is consistent with the implemented Trypsin/P rule.

8. **Medium: FASTA downloads are mutable and unverified.**  
   `get_fasta.sh` queries live UniProt and records neither release/date nor checksums. Re-running later can produce a different benchmark database. It also accepts an HTTP-successful but empty FASTA and still builds the combined file. Pin downloaded artifacts by checksum or archive them with UniProt release metadata; fail when counts are zero or expected suffix counts are implausible.

9. **Medium: DIA-NN installer can mark an unusable extraction successful.**  
   `get_diann.sh` touches `.ok` without checking the expected binary exists and is executable. An archive-layout change can create a permanently “successful” broken installation. Validate each version’s exact binary path, run a version/help probe, and only then create `.ok`. Add published archive checksums.

10. **Medium: `run.sh` performs insufficient input validation and unsafe config tokenization.**  
    It checks neither that each `.d` exists nor that at least one A and one B run was supplied. Also, `SETTINGS=$(...)` followed by unquoted `$SETTINGS` relies on shell word splitting and cannot safely represent quoted arguments. Parse the config into a Bash array, validate every run path, and check the output report exists and is nonempty after DIA-NN exits.

11. **Low: run filters can select unintended runs.**  
    Prefix `A1` also matches names beginning `A10`; substring slot matching can collide with longer identifiers. Support anchored regex or validate the complete selected run list and expected replicate count. A warning when filtering is disabled should preferably be a hard error for multi-run reports.

12. **Low: baseline rank sampling excludes curve index 0.**  
    `baseline_intensity` samples `[1, 10000)` then indexes the curve, so the highest-abundance first element is never selected. This is deterministic and does not affect A/B ratios, but it is not literally verbatim-equivalent to drawing among all 10,000 curve entries. Use an integer draw over `[0, len(curve))` if exact curve sampling is intended.
tokens used
14,676
**Findings**

1. **High: species assignment can become arbitrary during aggregation.**  
   `build_seed_from_report.py: groupby aggregation` selects `protein = s.iloc[0]`. If replicate rows report different protein groups, including cross-species assignments, the chosen species depends on row order. Imposed mode may therefore retain an ambiguous peptide and assign the wrong fraction, potentially differently in A and B. Aggregate all protein names for each `(sequence, charge)`, then call `assign_species` on the union; drop groups containing multiple species.

2. **High: partially missing charge silently deletes rows.**  
   `build_seed_from_report.py: has_charge/groupby` uses charge as a key whenever any charge is present. Pandas `groupby` drops rows whose charge is NaN, before the later optional-column logic can omit charge. Either require charge to be complete, use `dropna=False`, or explicitly separate/repair missing-charge rows. Also validate that charges are finite positive integers before casting to `Int64`.

3. **High: imposed mode is exact only on matched precursor ions, before independent seed normalization.**  
   For an ion present in both files with identical sequence, charge, and species assignment, the raw imposed ratio is exactly the composition ratio. TimSim’s independent median normalization changes it to:

   `observed B/A = nominal B/A × median(A seed intensities) / median(B seed intensities)`.

   A human-reference correction removes that common multiplier exactly in ideal arithmetic because every matched human ion has nominal ratio 1.0. The reasoning is therefore sound, but the claim should be qualified. Exact recovery fails or becomes empirical when:
   - A/B use different charge states or species assignments.
   - Normalization is estimated after missing values, stochastic simulation, censoring, saturation, or identification filtering.
   - There are too few shared human precursors.
   - “Human normalization” operates on protein summaries whose peptide sets differ.
   - Downstream processing applies nonlinear or run-specific normalization.

   Different A/B seed ID sets do not invalidate exactness for their intersection, but they can change both TimSim medians and the reference factor’s practical estimate.

4. **Medium: imposed mode unnecessarily depends on valid observed intensity.**  
   `build_seed_from_report.py: required-field filtering` drops rows with absent or nonpositive report intensity before imposed intensities are generated. This contradicts “discard found intensity”: an identified precursor with zero/missing quantity cannot enter imposed mode. In imposed mode, require identification fields but do not filter on reported intensity unless this is an intentional detectability criterion and is documented.

5. **Medium: sequence-level uniqueness hides charge-specific coverage and truth differences.**  
   `qc_fasta_coverage.py` and its reporting deduplicate only by `sequence`, while seeds and deterministic baselines are precursor-ion `(sequence, charge)` based. Peptide FASTA coverage is legitimately sequence-based, but report both unique peptides and affected seed ions; otherwise the apparent false-negative impact can be understated.

6. **Medium: QC digest default does not match the bundled DIA-NN search.**  
   `diann-pye1.cfg` specifies `--missed-cleavages 1`, while `qc_fasta_coverage.py` defaults to 2. This can label two-missed-cleavage seed peptides as searchable when DIA-NN will not search them. Change the default to 1 or pass search settings explicitly.

7. **Digest cleavage indexing itself is correct for Trypsin/P.**  
   `qc_fasta_coverage.py:digest` cuts after every K/R, including before proline, and the `j` bound correctly produces zero through `missed` missed cleavages. There is no off-by-one error. However, it is only an approximate DIA-NN search-space model: it does not validate amino-acid alphabet, terminal handling, precursor m/z/charge feasibility, modifications, or DIA-NN-specific peptide filtering. The configured `--cut K*,R*` is consistent with the implemented Trypsin/P rule.

8. **Medium: FASTA downloads are mutable and unverified.**  
   `get_fasta.sh` queries live UniProt and records neither release/date nor checksums. Re-running later can produce a different benchmark database. It also accepts an HTTP-successful but empty FASTA and still builds the combined file. Pin downloaded artifacts by checksum or archive them with UniProt release metadata; fail when counts are zero or expected suffix counts are implausible.

9. **Medium: DIA-NN installer can mark an unusable extraction successful.**  
   `get_diann.sh` touches `.ok` without checking the expected binary exists and is executable. An archive-layout change can create a permanently “successful” broken installation. Validate each version’s exact binary path, run a version/help probe, and only then create `.ok`. Add published archive checksums.

10. **Medium: `run.sh` performs insufficient input validation and unsafe config tokenization.**  
    It checks neither that each `.d` exists nor that at least one A and one B run was supplied. Also, `SETTINGS=$(...)` followed by unquoted `$SETTINGS` relies on shell word splitting and cannot safely represent quoted arguments. Parse the config into a Bash array, validate every run path, and check the output report exists and is nonempty after DIA-NN exits.

11. **Low: run filters can select unintended runs.**  
    Prefix `A1` also matches names beginning `A10`; substring slot matching can collide with longer identifiers. Support anchored regex or validate the complete selected run list and expected replicate count. A warning when filtering is disabled should preferably be a hard error for multi-run reports.

12. **Low: baseline rank sampling excludes curve index 0.**  
    `baseline_intensity` samples `[1, 10000)` then indexes the curve, so the highest-abundance first element is never selected. This is deterministic and does not affect A/B ratios, but it is not literally verbatim-equivalent to drawing among all 10,000 curve entries. Use an integer draw over `[0, len(curve))` if exact curve sampling is intended.
