"""Build a TimSim from_findings seed CSV from a real search-engine report.

Reads a DIA-NN or FragPipe report for ONE PYE1 sample and emits the 6-column
findings contract that TimSim's `from_findings` mode consumes:

    protein, sequence, intensity   (required)
    charge, rt, im                 (optional)

The A/B spike-in ratios are carried implicitly in the reported intensities and
organism membership — there is NO proteome_mix / dilution step. Run once per
sample (Sample A, Sample B).

Key contract details (see plasmabench/peptides.py and TimSim load_findings):
  - `sequence` must be UniMod-notated and parseable by the Rust PeptideSequence
    parser; rows that don't validate are dropped (with a count) so the Rust core
    never aborts.
  - Optional columns (rt/im) are written as real values or the column is OMITTED
    entirely — never blank strings. TimSim runs dropna() over every present
    column, so one all-empty optional column would delete every row.

Usage:
    python scripts/build_seed_from_report.py \
        --engine diann \
        --report /path/to/sampleA/report.tsv \
        --sample A \
        --out simulations/stage1/seeds/seed_sampleA.csv

    python scripts/build_seed_from_report.py \
        --engine fragpipe \
        --report /path/to/sampleB/psm.tsv \
        --sample B \
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
from plasmabench.species import assign_species  # noqa: E402

# Column name candidates per engine. Search-engine versions differ; the first
# present candidate wins, and we error clearly if none match (so a schema drift
# fails loudly instead of silently dropping data). Validate against Ute's
# actual exports — see TODO.md.
DIANN_COLS = {
    "protein": ["Protein.Group", "Protein.Ids", "Protein.Names"],
    "modseq": ["Modified.Sequence"],
    "charge": ["Precursor.Charge"],
    "intensity": ["Precursor.Quantity", "Precursor.Normalised", "Ms1.Area"],
    "rt_min": ["RT", "RT.Start"],          # DIA-NN reports RT in MINUTES
    "im": ["IM", "iIM", "Predicted.IM"],
}
FRAGPIPE_COLS = {
    "protein": ["Protein", "Protein ID"],
    "peptide": ["Peptide", "Peptide Sequence"],
    "mods": ["Assigned Modifications"],
    "charge": ["Charge"],
    "intensity": ["Intensity"],
    "rt_sec": ["Retention"],               # FragPipe psm.tsv Retention is SECONDS
    "im": ["Ion Mobility"],
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
    for required in ("protein", "modseq", "intensity"):
        if cols[required] is None:
            raise SystemExit(f"DIA-NN report missing a column for '{required}' "
                             f"(looked for {DIANN_COLS[required]})")
    out = pd.DataFrame({
        "protein": df[cols["protein"]].astype(str),
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
                    help="PYE1 sample label — used only for the run log / species audit")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    df = _read(args.report)
    out = load_diann(df) if args.engine == "diann" else load_fragpipe(df)
    n_raw = len(out)

    # Drop rows without the required fields (intensity > 0, real sequence).
    out = out[out["intensity"].notna() & (out["intensity"] > 0)]
    bad_seq = ~out["sequence"].map(is_valid_timsim_sequence)
    n_bad = int(bad_seq.sum())
    out = out[~bad_seq]

    # Aggregate to one row per (sequence, charge): sum intensity, median rt/im.
    has_charge = "charge" in out.columns and out["charge"].notna().any()
    group_keys = ["sequence", "charge"] if has_charge else ["sequence"]
    agg = {"intensity": ("intensity", "sum"),
           "protein": ("protein", lambda s: s.iloc[0])}
    if "rt" in out.columns:
        agg["rt"] = ("rt", "median")
    if "im" in out.columns:
        agg["im"] = ("im", "median")
    grouped = out.groupby(group_keys, as_index=False).agg(**agg)

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
    print(f"  engine={args.engine} sample={args.sample}")
    print(f"  {len(seed):,} rows  ({grouped['sequence'].nunique():,} unique sequences)"
          f"  from {n_raw:,} raw report rows")
    if n_bad:
        print(f"  dropped {n_bad:,} rows with sequences TimSim can't parse")
    if dropped_optional:
        print(f"  omitted optional column(s) {dropped_optional} (not fully populated)")
    print(f"  species split: {counts}")
    print("  expected PYE1 mass fractions for sample "
          f"{args.sample}: HUMAN .90 / "
          f"{'YEAST .02 / ECOLI .08' if args.sample == 'A' else 'YEAST .06 / ECOLI .04'}")


if __name__ == "__main__":
    main()
