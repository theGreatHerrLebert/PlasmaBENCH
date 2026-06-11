"""Calibrate Stage-2 YE seed intensities to the real plasma background (report-proxy).

The imposed YE seed has an arbitrary absolute scale. For superimpose-on-real-plasma we
want the simulated YE to sit at the level it occupies in the REAL PYE mix relative to the
human background it's laid on. Report-proxy approach (raw PYE .d unavailable):

  target_median[species] = plasma_human_median × (real_YE_median[species] / real_human_median)

where the YE/human relative scale comes from Ute's real PYE-mix DIA-NN report and the
plasma human median from the reference plasma's own DIA-NN search. One multiplicative
factor per species is applied to BOTH sample seeds, so the nominal A/B ratios (yeast 3.0,
E. coli 0.5) are preserved exactly — only the absolute level moves.

Caveat (documented): this matches report-quant medians (DIA-NN output domain), not raw
MS1/MS2 peak areas. Gold-standard calibration needs the raw PYE-mix .d.

Usage:
    python scripts/calibrate_stage2_seeds.py \
        --pye-report data/raw/dia/reports/diann/G1_nLC_tTOF_report.tsv \
        --plasma-report results/diann-2.5-plasma-hye/report.parquet \
        --seed-dir simulations/stage2/seeds --samples A B \
        --in-suffix _ye --out-suffix _ye_calib
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.species import assign_species

SPECIES = ["YEAST", "ECOLI"]


def _read(p: Path) -> pd.DataFrame:
    return pd.read_parquet(p) if p.suffix.lower() in (".parquet", ".pq") \
        else pd.read_csv(p, sep="\t", low_memory=False)


def real_relative_scale(pye_report: Path) -> dict:
    """median YE / median human precursor quantity in the real PYE mix, per species."""
    df = _read(pye_report)
    qcol = "Precursor.Quantity" if "Precursor.Quantity" in df else "Precursor.Normalised"
    df = df[(df["Q.Value"] <= 0.01) & (df[qcol] > 0)].copy()
    df["sp"] = df["Protein.Names"].map(assign_species)
    med = df.groupby(["Precursor.Id", "sp"])[qcol].median().reset_index()
    m = {s: med[med["sp"] == s][qcol].median() for s in ["HUMAN", *SPECIES]}
    return {s: m[s] / m["HUMAN"] for s in SPECIES}, m


def plasma_human_median(plasma_report: Path) -> float:
    df = _read(plasma_report)
    qcol = "Precursor.Normalised" if "Precursor.Normalised" in df else "Precursor.Quantity"
    df = df[(df["Q.Value"] <= 0.01) & (df[qcol] > 0)].copy()
    df["sp"] = df["Protein.Names"].map(assign_species)
    return float(df[df["sp"] == "HUMAN"][qcol].median())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pye-report", required=True, type=Path)
    ap.add_argument("--plasma-report", required=True, type=Path)
    ap.add_argument("--seed-dir", required=True, type=Path)
    ap.add_argument("--samples", nargs="+", default=["A", "B"])
    ap.add_argument("--in-suffix", default="_ye")
    ap.add_argument("--out-suffix", default="_ye_calib")
    args = ap.parse_args()

    rel, real_med = real_relative_scale(args.pye_report)
    H = plasma_human_median(args.plasma_report)
    target = {s: H * rel[s] for s in SPECIES}
    print(f"real PYE YE/human: " + ", ".join(f"{s} {rel[s]:.3f}" for s in SPECIES))
    print(f"plasma human median: {H:.3e}")
    print(f"target YE median:   " + ", ".join(f"{s} {target[s]:.3e}" for s in SPECIES))

    # pooled current YE median per species across all samples (one factor per species)
    seeds = {s_: pd.read_csv(args.seed_dir / f"seed_sample{s_}{args.in_suffix}.csv")
             for s_ in args.samples}
    for d in seeds.values():
        d["sp"] = d["protein"].map(assign_species)
    pooled = pd.concat(seeds.values(), ignore_index=True)
    k = {}
    for s in SPECIES:
        cur = pooled[pooled["sp"] == s]["intensity"].median()
        k[s] = target[s] / cur
        print(f"  {s}: current pooled median {cur:.3e}  → k={k[s]:.4f} ({np.log10(k[s]):+.2f} log10)")

    for s_, d in seeds.items():
        d = d.copy()
        d["intensity"] = d.apply(lambda r: r["intensity"] * k.get(r["sp"], 1.0), axis=1)
        d = d.drop(columns=["sp"])
        out = args.seed_dir / f"seed_sample{s_}{args.out_suffix}.csv"
        d.to_csv(out, index=False)
        # report calibrated medians + position vs human
        dd = pd.read_csv(out); dd["sp"] = dd["protein"].map(assign_species)
        pos = {s: dd[dd["sp"] == s]["intensity"].median() / H for s in SPECIES}
        print(f"wrote {out}  (calibrated YE/human median: "
              + ", ".join(f"{s} {pos[s]:.3f}" for s in SPECIES) + ")")


if __name__ == "__main__":
    main()
