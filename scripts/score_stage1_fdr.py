"""P5 — simulated-spike-in FDR / TPR for a Stage-1 run vs the TimSim blueprint.

Builds the run manifest from the SIM tree, loads both normalized tables, scores the
engine's found set against transmitted simulated truth per (sample, level), and writes
the FDR/recall panel. Reports the engine's self-claimed q-cut next to the *empirical*
FDR — the gap is the whole point.

Usage:
    python scripts/score_stage1_fdr.py \
        --report results/diann-2.5-AB-v3/report.parquet \
        --sim-dir simulations/stage1/dia \
        --experiment stage1-diann-2.5 \
        --out results/figures/stage1_diann25_fdr.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots.fdr import score_fdr, plot_fdr_tpr, blank_keysets


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="stage1")
    ap.add_argument("--human-scope", default="simulated",
                    choices=["simulated", "background_unknown"])
    ap.add_argument("--q-value-max", type=float, default=0.01)
    ap.add_argument("--blank-report", type=Path, default=None,
                    help="DIA-NN report of the BLANK noise source (same version/cfg). "
                         "Its IDs are subtracted from the FP side so blank contamination "
                         "doesn't inflate the empirical FDR.")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    manifest = build_manifest(args.report, args.sim_dir, args.experiment, args.human_scope)
    truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)
    obs = load_observations(args.report, manifest)

    blank_keys = None
    if args.blank_report is not None:
        blank_keys = blank_keysets(args.blank_report, q_value_max=args.q_value_max)
        print(f"blank subtraction: {{level: n}} = "
              f"{ {k: len(v) for k, v in blank_keys.items()} }")

    scores = score_fdr(truth, obs, q_value_max=args.q_value_max, blank_keys=blank_keys)
    print(f"human_scope={args.human_scope}  engine q-cut={args.q_value_max}\n")
    print(scores.to_string(index=False))

    plot_fdr_tpr(scores, args.out)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
