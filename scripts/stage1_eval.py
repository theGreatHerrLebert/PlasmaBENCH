"""Stage-1 evaluator — `make stage1-eval`. The EVIDENT replay command for the
`plasmabench-stage1-simulated-truth` claim; writes results/stage1/metrics.json.

Scores each DIA-NN report against the TimSim blueprint and emits species-resolved
ratio recovery + empirical FDR/recall. Follows the spec (memory stage1-eval-spec /
docs/codex-review-2026-06-08-e2e.md):
  1. Human-anchor the A/B ratios (build_ratio_table does this; both raw + anchored
     reported via the anchor value).
  2. Truth = the per-sample blueprint, intersected for ratio scoring (build_ratio_table
     keeps only features present in BOTH A and B).
  3. Singletons (present in one sample) are excluded from ratio scoring, counted in recall.
  4. CROSS-ENGINE: read differences, not absolute self-recall (seeding circularity).
  5. G-on-L hybrid: valid for engine-recovery scoring, NOT for absolute site/method claims.

Usage:
    python scripts/stage1_eval.py                 # defaults: 1.8 + 2.5 stage-1 reports
    python scripts/stage1_eval.py --report 1.8:results/diann-1.8-AB-v3/report.tsv ...
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots import ratios as R
from plasmabench.plots.fdr import score_fdr, blank_keysets
from plasmabench.species import EXPECTED_RATIO_B_OVER_A
import numpy as np

DEFAULT_REPORTS = ["1.8:results/diann-1.8-AB-v3/report.tsv",
                   "2.5:results/diann-2.5-AB-v3/report.parquet"]
RATIO_TOL_LOG2 = 0.4   # |observed - expected| A/B tolerance (imposed seeds, slight compression)
MIN_FEATURES = 500     # PASS requires a non-trivial scored set per level (spec: no tiny-set pass)
CAVEATS = [
    "Cross-engine: seeding from DIA-NN enriches truth for DIA-NN-findable peptides; read "
    "cross-engine DIFFERENCES, not absolute self-recall (spec req 4).",
    "G-on-L hybrid acquisition (site-G IDs, site-L window geometry): valid for engine-recovery "
    "scoring; do NOT make absolute site/acquisition-method sensitivity claims (spec req 5).",
    "Ratios are human-anchored (one run-pair scale from HUMAN precursors); the raw blueprint "
    "carries a uniform offset from per-seed median normalization (spec req 1).",
]


def score_one(label: str, report: Path, sim_dir: Path, blank_keys, q: float) -> dict:
    manifest = build_manifest(report, sim_dir, f"stage1-{label}", "simulated")
    truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)
    obs = load_observations(report, manifest)

    # spec req 2/3: ratio recovery is scored ONLY on observed precursors that are in the
    # TRANSMITTED blueprint (excludes false positives from the ratio oracle; singletons are
    # already dropped by build_ratio_table's both-sample requirement).
    # Per-SAMPLE key (sample, seq, charge): an ion transmitted only in A must not credit a
    # false B observation (a global seq+charge set would). Both members of an A/B pair must
    # each be a transmitted-blueprint precursor in their own sample.
    tt = truth[truth["transmitted"]]
    truth_keys = set(zip(tt["sample"], tt["sequence_modified"], tt["charge"]))
    okeys = list(zip(obs["sample"], obs["sequence_modified"], obs["charge"]))
    obs_truth = obs[[k in truth_keys for k in okeys]].copy()

    # --- ratio recovery (human-anchored; blueprint-intersected; both-sample only) ---
    ratio = {}
    for level, pq in (("ion", "maxlfq"), ("protein", "raw_sum")):
        wide, anchor = R.build_ratio_table(obs_truth, level=level, q_value_max=q,
                                           human_anchor=True, protein_quant=pq)
        summ = R.ratio_summary(wide, R.A_OVER_B).set_index("species")
        per_sp = {}; worst = 0.0
        for sp in ("YEAST", "ECOLI"):
            obs_l2 = float(summ.loc[sp, "observed_log2"]); exp_l2 = float(summ.loc[sp, "expected_log2"])
            err = abs(obs_l2 - exp_l2); worst = max(worst, err)
            per_sp[sp] = {"observed_log2_AB": round(obs_l2, 3),
                          "expected_log2_AB": round(exp_l2, 3), "abs_error_log2": round(err, 3)}
        ratio[level] = {"n_features": int(len(wide)), "human_anchor_log2": round(float(anchor), 3),
                        "species": per_sp, "worst_abs_error_log2": round(worst, 3),
                        "within_tolerance": bool(worst <= RATIO_TOL_LOG2 and len(wide) >= MIN_FEATURES)}

    # --- empirical FDR / recall vs blueprint (blank-subtracted); POOL TP/FP/FN across
    #     samples (one confusion matrix), not a macro-average of per-sample percentages ---
    sc = score_fdr(truth, obs, q_value_max=q, blank_keys=blank_keys)
    fdr_recall = {}
    for level in ("ion", "peptide", "protein"):
        sub = sc[sc["level"] == level]
        if level == "protein":
            sub = sub[sub["mode"] == "group-credit:any"]
        if len(sub):
            tp, fp = float(sub["TP"].sum()), float(sub["FP"].sum())
            fdr = round(100.0 * fp / (fp + tp), 3) if (fp + tp) else 0.0
            if level == "protein":
                # protein TP/FP are GROUP counts (group-credit); recall must use protein counts,
                # not TP/(TP+FN) (FN is a truth-protein count — incompatible unit).
                dp, ntp = float(sub["detected_proteins"].sum()), float(sub["n_truth_proteins"].sum())
                recall = round(100.0 * dp / ntp, 3) if ntp else 0.0
            else:
                fn = float(sub["FN"].sum())
                recall = round(100.0 * tp / (tp + fn), 3) if (tp + fn) else 0.0
            fdr_recall[level] = {"fdr_pct": fdr, "recall_pct": recall, "TP": int(tp), "FP": int(fp)}
    return {"ratio_recovery": ratio, "fdr_recall": fdr_recall,
            "n_truth_precursors": int((truth["transmitted"]).sum())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="append", default=None,
                    help="label:path (repeatable). Default: 1.8 + 2.5 stage-1 reports.")
    ap.add_argument("--sim-dir", type=Path, default=Path("simulations/stage1/dia"))
    ap.add_argument("--blank-report", type=Path, default=Path("results/diann-2.5-blank/report.parquet"))
    ap.add_argument("--q-value-max", type=float, default=0.01)
    ap.add_argument("--out", type=Path, default=Path("results/stage1/metrics.json"))
    args = ap.parse_args()
    reports = args.report or DEFAULT_REPORTS

    missing = [r.split(":", 1)[1] for r in reports if not Path(r.split(":", 1)[1]).exists()]
    if missing or not args.sim_dir.exists():
        print(f"stage1-eval: required inputs missing (reports {missing}, sim-dir {args.sim_dir}). "
              f"Run the Stage-1 SIM + DIA-NN searches first (see CLAUDE.md).", file=sys.stderr)
        return 65

    blank_keys = None
    if args.blank_report.exists():
        blank_keys = blank_keysets(args.blank_report, q_value_max=args.q_value_max)
    else:
        print(f"WARNING: blank report {args.blank_report} not found — FDR is NOT blank-subtracted "
              f"(blank-leaked IDs may inflate it). Pass --blank-report to subtract.", file=sys.stderr)

    engines = {}
    for spec in reports:
        label, path = spec.split(":", 1)
        print(f"scoring DIA-NN {label} ...", file=sys.stderr)
        engines[label] = score_one(label, Path(path), args.sim_dir, blank_keys, args.q_value_max)

    # cross-engine DELTAS — the unbiased read (spec req 4): absolute self-recall is biased by
    # seeding circularity, so report engine-to-engine differences, not absolutes.
    labels = list(engines)
    cross = {}
    if len(labels) == 2:
        a, b = labels
        def delta(path_fn):
            return round(path_fn(engines[a]) - path_fn(engines[b]), 3)
        cross = {
            "engines": [a, b],
            "ratio_ion_worst_abs_error_log2_delta": delta(
                lambda e: e["ratio_recovery"]["ion"]["worst_abs_error_log2"]),
            "ion_fdr_pct_delta": delta(lambda e: e["fdr_recall"].get("ion", {}).get("fdr_pct", 0.0)),
            "ion_recall_pct_delta": delta(lambda e: e["fdr_recall"].get("ion", {}).get("recall_pct", 0.0)),
            "protein_fdr_pct_delta": delta(lambda e: e["fdr_recall"].get("protein", {}).get("fdr_pct", 0.0)),
            "note": f"deltas = {a} minus {b}; read differences, not absolute self-recall.",
        }

    overall_pass = all(e["ratio_recovery"]["ion"]["within_tolerance"] and
                       e["ratio_recovery"]["protein"]["within_tolerance"] for e in engines.values())
    metrics = {
        "schema": "plasmabench.stage1.metrics/0.1",
        "decision_rule": f"PASS iff every engine's yeast & E.coli A/B (human-anchored, ion AND "
                         f"protein) recover nominal within {RATIO_TOL_LOG2} log2, each level "
                         f"scoring >= {MIN_FEATURES} blueprint-intersected features.",
        "result": "PASS" if overall_pass else "FAIL",
        "sim_dir": str(args.sim_dir),
        "q_value_max": args.q_value_max,
        "blank_subtracted": blank_keys is not None,
        "expected_ratio_B_over_A": EXPECTED_RATIO_B_OVER_A,
        "caveats": CAVEATS,
        "engines": engines,
        "cross_engine": cross,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(metrics, indent=2))

    # human-readable summary
    print(f"\n=== Stage-1 evaluation  ({metrics['result']}) ===")
    for la, e in engines.items():
        ri = e["ratio_recovery"]["ion"]["species"]; rp = e["ratio_recovery"]["protein"]["species"]
        fr = e["fdr_recall"]
        print(f"  DIA-NN {la}: ion A/B ecoli {ri['ECOLI']['observed_log2_AB']:+.2f} / yeast "
              f"{ri['YEAST']['observed_log2_AB']:+.2f}  | ion FDR {fr.get('ion',{}).get('fdr_pct','-')}% "
              f"recall {fr.get('ion',{}).get('recall_pct','-')}%  | protein FDR "
              f"{fr.get('protein',{}).get('fdr_pct','-')}%")
    print(f"wrote {args.out}")
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
