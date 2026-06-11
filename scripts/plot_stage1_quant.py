"""P6 quant-correlation for a Stage-1 run: simulated truth vs DIA-NN quantity.

Builds the run manifest from the SIM output tree (each run's blueprint lives at
``<sim-dir>/**/<run_id>/synthetic_data.db``, where ``run_id`` == DIA-NN ``Run``),
loads the two normalized tables, and writes the precursor-level correlation panel.

Usage:
    python scripts/plot_stage1_quant.py \
        --report results/diann-2.5-AB-v3/report.parquet \
        --sim-dir simulations/stage1/dia \
        --experiment stage1-diann-2.5 \
        --out results/figures/stage1_diann25_quant.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots.quant import build_quant_table, correlations, plot_quant_correlation


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="stage1")
    ap.add_argument("--human-scope", default="simulated",
                    choices=["simulated", "background_unknown"])
    ap.add_argument("--q-value-max", type=float, default=0.01)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    manifest = build_manifest(args.report, args.sim_dir, args.experiment, args.human_scope)
    print(f"manifest: {len(manifest.runs)} runs")
    for r in manifest.runs:
        print(f"  {r.run_id}  sample={r.sample}  -> {r.blueprint_db}")

    truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)
    obs = load_observations(args.report, manifest)
    print(f"truth: {len(truth):,} precursors ({int(truth['transmitted'].sum()):,} transmitted)  "
          f"observations: {len(obs):,}")

    merged = build_quant_table(truth, obs, q_value_max=args.q_value_max)
    corrs = correlations(merged)
    print(f"\nmatched {len(merged):,} precursors (transmitted truth ∩ q<{args.q_value_max} obs)")
    print(corrs.to_string(index=False))

    plot_quant_correlation(merged, corrs, args.out)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
