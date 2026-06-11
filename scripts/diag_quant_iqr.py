"""Diagnostic: A/B ratio IQR per engine, FULL vs SHARED-precursor intersection.

Settles the coverage-confound question Codex raised: is DIA-NN 2.5's wider A/B
ratio IQR a property of its quant on the SAME precursors, or just an artifact of
it reporting more (marginal) precursors? Compares two engines' ion-level ratio
tables on (a) each engine's full precursor set and (b) the precursors they share.

Usage:
    python scripts/diag_quant_iqr.py \
        --engine 1.8 results/diann-1.8-AB-v3/report.parquet \
        --engine 2.5 results/diann-2.5-AB-v3/report.parquet \
        --sim-dir simulations/stage1/dia --experiment stage1 --human-scope simulated
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots.ratios import build_ratio_table, SPECIES_ORDER


def iqr(x: pd.Series) -> float:
    x = x.dropna()
    if len(x) < 4:
        return float("nan")
    return float(np.subtract(*np.percentile(x, [75, 25])))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", nargs=2, action="append", metavar=("LABEL", "REPORT"),
                    required=True, help="repeat twice: --engine 1.8 path --engine 2.5 path")
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="diag")
    ap.add_argument("--human-scope", default="simulated",
                    choices=["simulated", "background_unknown"])
    ap.add_argument("--level", default="ion")
    ap.add_argument("--q-value-max", type=float, default=0.01)
    args = ap.parse_args()

    tables: dict[str, pd.DataFrame] = {}
    for label, report in args.engine:
        report = Path(report)
        manifest = build_manifest(report, args.sim_dir, args.experiment, args.human_scope)
        obs = load_observations(report, manifest)
        wide, anchor = build_ratio_table(obs, level=args.level, q_value_max=args.q_value_max)
        wide = wide.set_index("feature")
        tables[label] = wide
        print(f"[{label}] {report}  anchor={anchor:+.3f}  features={len(wide):,}")

    labels = [lab for lab, _ in args.engine]
    a, b = tables[labels[0]], tables[labels[1]]

    print(f"\n{'species':8} {'engine':6} {'n_full':>8} {'IQR_full':>9} "
          f"{'n_shared':>9} {'IQR_shared':>11}")
    for sp in [s for s in SPECIES_ORDER if s != "HUMAN"]:
        ai = set(a.index[a["species"] == sp])
        bi = set(b.index[b["species"] == sp])
        shared = ai & bi
        for lab, tbl in zip(labels, (a, b)):
            full = tbl.loc[tbl["species"] == sp, "log2_ba"]
            sh = tbl.loc[tbl.index.isin(shared), "log2_ba"]
            print(f"{sp:8} {lab:6} {len(full):>8} {iqr(full):>9.3f} "
                  f"{len(sh):>9} {iqr(sh):>11.3f}")
        print(f"{'':8} shared set: {len(shared):,}  "
              f"({labels[0]} only {len(ai - bi):,}, {labels[1]} only {len(bi - ai):,})")


if __name__ == "__main__":
    main()
