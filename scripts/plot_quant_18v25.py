"""Report figures for the DIA-NN 1.8 vs 2.5 quant finding (report §8).

Panel A: A/B ratio precision (IQR of log2 fold change) by engine × stage × species,
         using the default Precursor.Normalised quantity (what a user gets out-of-the-box).
Panel B: mechanism — IQR by quantity column for Stage 2 (1.8 vs 2.5). Shows 2.5's spread
         is the cross-run normalization: raw Precursor.Quantity / Ms1.Area recover 1.8 levels.

Reuses the loader + ratio table; computes everything (no hardcoded numbers).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_observations
from plasmabench.plots.ratios import build_ratio_table, SPECIES_COLOR

FIGDIR = Path("results/figures")

STAGES = {
    "Stage 1\n(blank)": dict(sim="simulations/stage1/dia", exp="stage1", scope="simulated",
                             reps={"1.8": "results/diann-1.8-AB-v3/report.tsv",
                                   "2.5": "results/diann-2.5-AB-v3/report.parquet",
                                   "2.6": "results/diann-2.6-AB/report.parquet"}),
    "Stage 2\n(real plasma)": dict(sim="simulations/stage2/dia-fix", exp="stage2",
                                   scope="background_unknown",
                                   reps={"1.8": "results/diann-1.8-s2fix-r080/report.tsv",
                                         "2.5": "results/diann-2.5-s2fix-r080/report.parquet",
                                         "2.6": "results/diann-2.6-s2fix-r080/report.parquet"}),
}
QUANT_COLS = ["Precursor.Normalised", "Precursor.Quantity", "Ms1.Normalised", "Ms1.Area"]
SPECIES = ["YEAST", "ECOLI"]
ENGINES = ["1.8", "2.5", "2.6"]
COLORS = {"1.8": "#4C72B0", "2.5": "#C44E52", "2.6": "#DD8452"}


def iqr(x):
    x = pd.Series(x).dropna()
    return float(np.subtract(*np.percentile(x, [75, 25]))) if len(x) >= 4 else np.nan


def species_iqr(report, sim, exp, scope, quant_col=None):
    manifest = build_manifest(Path(report), Path(sim), exp, scope)
    obs = load_observations(Path(report), manifest, quant_col=quant_col)
    wide, _ = build_ratio_table(obs, level="ion", q_value_max=0.01, human_anchor=True)
    return {sp: iqr(wide.loc[wide["species"] == sp, "log2_ba"]) for sp in SPECIES}


def panel_a():
    rows = []
    for stage, cfg in STAGES.items():
        for eng, rep in cfg["reps"].items():
            d = species_iqr(rep, cfg["sim"], cfg["exp"], cfg["scope"],
                            quant_col="Precursor.Quantity")
            for sp in SPECIES:
                rows.append((stage, eng, sp, d[sp]))
    df = pd.DataFrame(rows, columns=["stage", "engine", "species", "iqr"])

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    stages = list(STAGES)
    for ax, stage in zip(axes, stages):
        sub = df[df["stage"] == stage]
        x = np.arange(len(SPECIES)); w = 0.27
        for i, eng in enumerate(ENGINES):
            vals = [sub[(sub.engine == eng) & (sub.species == sp)]["iqr"].iloc[0] for sp in SPECIES]
            bars = ax.bar(x + (i - (len(ENGINES) - 1) / 2) * w, vals, w,
                          label=f"DIA-NN {eng}", color=COLORS[eng])
            for b, v in zip(bars, vals):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                        ha="center", va="bottom", fontsize=7)
        ax.set_title(stage.replace("\n", " "), fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels([s.title() for s in SPECIES])
        ax.axhline(0, color="0.6", lw=0.6)
    axes[0].set_ylabel("A/B ratio IQR  (log2; lower = tighter)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("A/B ratio precision — raw Precursor.Quantity (benchmark standard): 1.8 ≈ 2.5 ≈ 2.6",
                 fontsize=11)
    fig.tight_layout()
    out = FIGDIR / "quant_18v25_iqr.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def panel_b():
    cfg = STAGES["Stage 2\n(real plasma)"]
    rows = []
    for eng, rep in cfg["reps"].items():
        avail = (pd.read_parquet(rep).columns if str(rep).endswith(".parquet")
                 else pd.read_csv(rep, sep="\t", nrows=1).columns)
        for qc in [c for c in QUANT_COLS if c in set(avail)]:
            d = species_iqr(rep, cfg["sim"], cfg["exp"], cfg["scope"], quant_col=qc)
            for sp in SPECIES:
                rows.append((eng, qc, sp, d[sp]))
    df = pd.DataFrame(rows, columns=["engine", "quant", "species", "iqr"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    cols = [c for c in QUANT_COLS if c in set(df["quant"])]
    x = np.arange(len(cols)); w = 0.27
    for ax, sp in zip(axes, SPECIES):
        for i, eng in enumerate(ENGINES):
            vals = [df[(df.engine == eng) & (df.quant == qc) & (df.species == sp)]["iqr"]
                    for qc in cols]
            vals = [v.iloc[0] if len(v) else np.nan for v in vals]
            off = (i - (len(ENGINES) - 1) / 2) * w
            ax.bar(x + off, vals, w, label=f"DIA-NN {eng}", color=COLORS[eng])
            for xi, v in zip(x + off, vals):
                if not np.isnan(v):
                    ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5)
        ax.set_title(sp.title(), fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace(".", ".\n") for c in cols], fontsize=8)
    axes[0].set_ylabel("A/B ratio IQR  (log2; lower = tighter)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Stage 2: the cross-run NORMALIZATION spread (2.5 AND 2.6, not 1.8) "
                 "— raw Quantity / Ms1.Area recover 1.8 levels", fontsize=10)
    fig.tight_layout()
    out = FIGDIR / "quant_18v25_normalization.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    FIGDIR.mkdir(parents=True, exist_ok=True)
    panel_a()
    panel_b()
