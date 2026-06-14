"""Capstone: real-vs-SIM A/B ratio vs abundance (the simulator low-input over-separation gap).

Overlays the ratio-vs-abundance curve, per species, for:
  - REAL G PYE1 (experimental), as the MEAN of the 6 single-A-vs-single-B injection pairs
    — matched to SIM's 1-vs-1 structure (1 A run, 1 B run). This matters: averaging all 6
    real replicates per condition (min_replicates=1, linear_mean) artificially stabilises
    low-abundance ratios and manufactures apparent "compression"; the matched 1-vs-1 curve
    removes that confound (Codex review, 2026-06-14).
  - SIM Stage 2 (YE spike-in on the SAME real plasma batch) — already 1-vs-1.
  - SIM Stage 1 (blank background) — already 1-vs-1.
All RAW Precursor.Quantity, ion level, human-anchored, DIA-NN 1.8.

x-axis = abundance OCTILE RANK (1=lowest…8=highest); curves compare each dataset's own
low-abundance tail. CAVEAT: rank is not a matched absolute input across datasets (different
intensity scales), so the claim is the within-dataset SHAPE: SIM bends AWAY from nominal
(over-separates) at low rank; REAL stays flat. Mechanism (detector vs imposed-intensity) is
not isolated here.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import (build_manifest, load_observations,
                                       load_observations_real)
from plasmabench.plots import ratios as R

FIG = Path("results/figures/real_vs_sim_compression.png")
NBINS = 8
SPECIES = ["ECOLI", "YEAST", "HUMAN"]
EXP = {sp: R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA[sp]) for sp in SPECIES}


def octile_curve(obs):
    w, _ = R.build_ratio_table(obs, level="ion", q_value_max=0.01, human_anchor=True)
    b = R.bias_by_abundance(w, R.A_OVER_B, n_bins=NBINS)
    return {sp: b[b.species == sp].sort_values("x")["median"].to_numpy() for sp in SPECIES}


def real_matched_curve():
    """REAL as the mean over the 6 single-A vs single-B injection pairs (1-vs-1, matches SIM)."""
    obs = load_observations_real("results/diann-1.8-realG/report.tsv", quant_col="Precursor.Quantity")
    A = sorted(obs.loc[obs["sample"] == "A", "run_id"].unique())
    B = sorted(obs.loc[obs["sample"] == "B", "run_id"].unique())
    pairs = [octile_curve(obs[obs["run_id"].isin([a, b])]) for a, b in zip(A, B)]
    return {sp: np.mean([p[sp] for p in pairs], axis=0) for sp in SPECIES}, len(pairs)


def sim_curve(report, sim_dir, exp, scope):
    return octile_curve(load_observations(report, build_manifest(report, sim_dir, exp, scope)))


def main():
    FIG.parent.mkdir(parents=True, exist_ok=True)
    real, npairs = real_matched_curve()
    data = [
        (f"REAL G PYE1 (mean of {npairs} 1-vs-1 pairs)", "#222222", "o", "-", real),
        ("SIM Stage 2 (real plasma)", "#DD8452", "s", "--",
         sim_curve("results/diann-1.8-s2fix-r080/report.tsv", "simulations/stage2/dia-fix", "s2", "background_unknown")),
        ("SIM Stage 1 (blank)", "#4C72B0", "^", ":",
         sim_curve("results/diann-1.8-AB-v3/report.tsv", "simulations/stage1/dia", "s1", "simulated")),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, sp in zip(axes, SPECIES):
        for name, c, mk, ls, cur in data:
            y = cur[sp]
            ax.plot(np.arange(1, len(y) + 1), y, ls, color=c, marker=mk, ms=5, lw=1.8, label=name)
        ax.axhline(EXP[sp], color="0.5", lw=1.0, ls=(0, (4, 3)))
        ax.text(0.98, 0.93 if sp == "ECOLI" else 0.06, f"nominal {EXP[sp]:+.2f}",
                transform=ax.transAxes, ha="right", fontsize=8, color="0.4")
        ax.set_title(sp.title() + ("  (control)" if sp == "HUMAN" else ""), fontsize=11)
        ax.set_xlabel("abundance octile  (1 = lowest → 8 = highest)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("A/B ratio  (log2, human-anchored)")
    axes[0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Low-input over-separation is a SIMULATOR artifact: SIM ratios blow outward at "
                 "low abundance; matched real data stays flat (DIA-NN 1.8, raw quant)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"wrote {FIG}")
    for sp in SPECIES:
        print(f"\n{sp} (octile 1→8):")
        for name, *_, cur in data:
            print(f"  {name:34} " + " ".join(f"{v:+.2f}" for v in cur[sp]))


if __name__ == "__main__":
    main()
