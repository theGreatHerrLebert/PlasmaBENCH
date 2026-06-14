"""Capstone: real-vs-SIM A/B ratio compression vs abundance (the simulator-gap figure).

Overlays the ratio-bias-vs-abundance curve, per species, for:
  - REAL G PYE1 (Ute's experimental runs, our unified search)
  - SIM Stage 2 (YE spike-in superimposed on the SAME real plasma batch)
  - SIM Stage 1 (blank background)
all on RAW Precursor.Quantity, ion level, human-anchored, DIA-NN 1.8 (cleanest quant).

x-axis = abundance OCTILE RANK (1=lowest … 8=highest) so real and SIM align despite
different absolute intensity scales. The finding: REAL bends toward 0 (compresses) at low
rank; SIM bends AWAY from 0 (over-separates) at low rank — opposite low-input behaviour.
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
EXP = {"ECOLI": R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA["ECOLI"]),
       "YEAST": R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA["YEAST"]),
       "HUMAN": 0.0}
DATASETS = [
    ("REAL G PYE1", "#222222", "o", "-",
     lambda: load_observations_real("results/diann-1.8-realG/report.tsv",
                                    quant_col="Precursor.Quantity")),
    ("SIM Stage 2 (real plasma)", "#DD8452", "s", "--",
     lambda: load_observations("results/diann-1.8-s2fix-r080/report.tsv",
                               build_manifest("results/diann-1.8-s2fix-r080/report.tsv",
                                              "simulations/stage2/dia-fix", "s2",
                                              "background_unknown"))),
    ("SIM Stage 1 (blank)", "#4C72B0", "^", ":",
     lambda: load_observations("results/diann-1.8-AB-v3/report.tsv",
                               build_manifest("results/diann-1.8-AB-v3/report.tsv",
                                              "simulations/stage1/dia", "s1", "simulated"))),
]


def curve(obs):
    w, _ = R.build_ratio_table(obs, level="ion", q_value_max=0.01, human_anchor=True)
    b = R.bias_by_abundance(w, R.A_OVER_B, n_bins=NBINS)
    out = {}
    for sp in SPECIES:
        d = b[b.species == sp].sort_values("x")
        out[sp] = d["median"].to_numpy()
    return out


def main():
    FIG.parent.mkdir(parents=True, exist_ok=True)
    data = [(name, c, mk, ls, curve(load())) for name, c, mk, ls, load in DATASETS]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, sp in zip(axes, SPECIES):
        for name, c, mk, ls, cur in data:
            y = cur[sp]
            ax.plot(np.arange(1, len(y) + 1), y, ls, color=c, marker=mk, ms=5, lw=1.8,
                    label=name)
        ax.axhline(EXP[sp], color="0.5", lw=1.0, ls=(0, (4, 3)))
        ax.text(0.98, 0.03 if sp != "ECOLI" else 0.93, f"nominal {EXP[sp]:+.2f}",
                transform=ax.transAxes, ha="right", fontsize=8, color="0.4")
        ax.set_title(sp.title() + ("  (control)" if sp == "HUMAN" else ""), fontsize=11)
        ax.set_xlabel("abundance octile  (1 = lowest → 8 = highest)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("A/B ratio  (log2, human-anchored)")
    axes[0].legend(frameon=False, fontsize=8.5, loc="best")
    fig.suptitle("Real vs SIM: spike-in ratio compression at low abundance — "
                 "REAL compresses toward 1:1, SIM over-separates (DIA-NN 1.8, raw quant)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"wrote {FIG}")
    # also dump the numbers
    for sp in SPECIES:
        print(f"\n{sp} (octile 1→8):")
        for name, *_ , cur in data:
            print(f"  {name:28} " + " ".join(f"{v:+.2f}" for v in cur[sp]))


if __name__ == "__main__":
    main()
