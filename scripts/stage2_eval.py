"""Stage-2 evaluator — `make stage2-eval`. The EVIDENT replay for the
`plasmabench-stage2-superimpose-truth` claim; writes results/stage2/metrics.json
AND the cited ratio panel.

VERIFIES the claim's tolerance: the corrected-sim BLUEPRINT B/A is within 5% of
nominal (yeast 3.0, E. coli 0.5) on BOTH reference-plasma backgrounds (r080, r081)
— the sim-QC gate after the from_findings fix. Also reports the observed (DIA-NN 1.8)
ratio recovery at ion + protein level (raw quantity, human-anchored) as descriptive.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots import ratios as R

SIM_DIR = Path("simulations/stage2/dia-fix")
BACKGROUNDS = {"r080": "results/diann-1.8-s2fix-r080/report.tsv",
               "r081": "results/diann-1.8-s2fix-r081/report.tsv"}
NOMINAL_BA = {"YEAST": 3.0, "ECOLI": 0.5}
SIMQC_REL_TOL = 0.05   # blueprint B/A within 5% of nominal — the claim's gate
FIG = Path("results/figures/stage2fix_18_ratios_protein.png")
OUT = Path("results/stage2/metrics.json")


def blueprint_ba(truth: pd.DataFrame) -> dict:
    t = truth[truth["truth_intensity"] > 0]
    s = t.groupby(["species", "sample"])["truth_intensity"].sum().unstack()
    out = {}
    for sp in ("YEAST", "ECOLI"):
        if sp in s.index and {"A", "B"} <= set(s.columns):
            out[sp] = float(s.loc[sp, "B"] / s.loc[sp, "A"])
    return out


def main() -> int:
    missing = [p for p in BACKGROUNDS.values() if not Path(p).exists()] + \
              ([str(SIM_DIR)] if not SIM_DIR.exists() else [])
    if missing:
        print(f"stage2-eval: required inputs missing: {missing}. Run the Stage-2 SIM + DIA-NN "
              f"searches first (see CLAUDE.md).", file=sys.stderr)
        return 65

    backgrounds = {}
    simqc_pass = True
    for bg, report in BACKGROUNDS.items():
        m = build_manifest(Path(report), SIM_DIR, f"stage2-{bg}", "background_unknown")
        truth = pd.concat([load_truth(r, m) for r in m.runs], ignore_index=True)
        ba = blueprint_ba(truth)
        simqc = {}
        for sp, val in ba.items():
            rel = abs(val - NOMINAL_BA[sp]) / NOMINAL_BA[sp]
            ok = rel <= SIMQC_REL_TOL
            simqc_pass = simqc_pass and ok
            simqc[sp] = {"blueprint_BA": round(val, 3), "nominal_BA": NOMINAL_BA[sp],
                         "rel_error": round(rel, 3), "within_5pct": bool(ok)}
        # observed (DIA-NN 1.8) recovery, raw quant, human-anchored (descriptive)
        obs = load_observations(Path(report), m)
        observed = {}
        for level, pq in (("ion", "maxlfq"), ("protein", "raw_sum")):
            wide, anchor = R.build_ratio_table(obs, level=level, protein_quant=pq, human_anchor=True)
            summ = R.ratio_summary(wide, R.A_OVER_B).set_index("species")
            observed[level] = {sp: round(float(summ.loc[sp, "observed_log2"]), 3)
                               for sp in ("YEAST", "ECOLI")}
            if bg == "r080" and level == "protein":
                R.plot_ratio_panels(wide, FIG, conv=R.A_OVER_B, level=level, anchor=anchor)
        backgrounds[bg] = {"sim_qc": simqc, "observed_recovery_AB": observed}

    metrics = {
        "schema": "plasmabench.stage2.metrics/0.1",
        "decision_rule": f"PASS iff blueprint B/A is within {int(SIMQC_REL_TOL*100)}% of nominal "
                         f"(yeast 3.0, E. coli 0.5) on BOTH backgrounds (r080, r081).",
        "result": "PASS" if simqc_pass else "FAIL",
        "nominal_BA": NOMINAL_BA,
        "caveats": [
            "Human is real plasma background (no blueprint) — scored over the simulated spike-in only.",
            "Superimposition does not reproduce ion suppression / fill-time / detector saturation.",
        ],
        "backgrounds": backgrounds,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(metrics, indent=2))

    print(f"\n=== Stage-2 evaluation  ({metrics['result']}) ===")
    for bg, b in backgrounds.items():
        q = b["sim_qc"]
        print(f"  {bg}: blueprint B/A yeast {q.get('YEAST',{}).get('blueprint_BA','-')} "
              f"ecoli {q.get('ECOLI',{}).get('blueprint_BA','-')}  (nominal 3.0 / 0.5)")
    print(f"wrote {OUT} and {FIG}")
    return 0 if simqc_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
