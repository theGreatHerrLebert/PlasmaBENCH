"""LFQbench-style three-species ratio figure for a PlasmaBENCH stage-1 run.

Reads a DIA-NN report (parquet or tsv), computes per-precursor log2(B/A) from the
normalised quantities, assigns species, and draws the classic mixed-proteome
benchmark figure: (left) per-species log2-ratio distributions vs the expected
PYE1 ratios, (right) the LFQbench "MA" scatter of log2(B/A) vs abundance.

Ratios are human-anchored (the designed 1:1 background): all log2 ratios are
shifted so the human median sits at 0, which is standard for mixed-proteome
benchmarks and matches the stage-1 scorer's human-anchor normalisation.

Usage:
    python scripts/plot_stage1_ratios.py \
        --report results/diann-2.5-AB-v3/report.parquet \
        --out results/figures/stage1_diann25_ratios.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.species import assign_species, EXPECTED_RATIO_B_OVER_A  # noqa: E402

# --- timsim-bench house style ----------------------------------------------
FONT_AXIS_LABEL, FONT_TICK_LABEL, FONT_TITLE, FONT_ANNOTATION, FONT_LEGEND = 14, 12, 14, 14, 12
SPECIES_COLOR = {"HUMAN": "#4C72B0", "YEAST": "#DD8452", "ECOLI": "#55A868"}
SPECIES_ORDER = ["ECOLI", "HUMAN", "YEAST"]
SPECIES_LABEL = {"HUMAN": "Human (plasma)", "YEAST": "Yeast", "ECOLI": "E. coli"}


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() in (".parquet", ".pq") \
        else pd.read_csv(path, sep="\t", low_memory=False)


def build_ratios(df: pd.DataFrame) -> pd.DataFrame:
    qcol = "Precursor.Normalised" if "Precursor.Normalised" in df.columns else "Precursor.Quantity"
    pname = "Protein.Names" if "Protein.Names" in df.columns else "Protein.Group"
    df = df.copy()
    df["sample"] = np.where(df["Run"].str.contains("sampleA"), "A",
                            np.where(df["Run"].str.contains("sampleB"), "B", None))
    df = df[df["sample"].notna() & (df[qcol] > 0)]
    df["species"] = df[pname].map(assign_species)
    wide = df.pivot_table(index=["Precursor.Id", "species"], columns="sample",
                          values=qcol, aggfunc="median")
    wide = wide.dropna(subset=["A", "B"]).reset_index()
    wide = wide[wide["species"].isin(SPECIES_ORDER)]
    wide["log2BA"] = np.log2(wide["B"] / wide["A"])
    wide["log2int"] = 0.5 * np.log2(wide["A"] * wide["B"])
    # human-anchor: shift so the 1:1 background centres at 0
    anchor = wide.loc[wide["species"] == "HUMAN", "log2BA"].median()
    wide["log2BA"] -= anchor
    return wide, anchor


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--title", default="PlasmaBENCH Stage 1 — DIA-NN 2.5 on simulated PYE1")
    args = ap.parse_args()

    w, anchor = build_ratios(_read(args.report))
    exp_log2 = {s: np.log2(EXPECTED_RATIO_B_OVER_A[s]) for s in SPECIES_ORDER}

    mpl.rcParams.update({"axes.labelsize": FONT_AXIS_LABEL, "xtick.labelsize": FONT_TICK_LABEL,
                         "ytick.labelsize": FONT_TICK_LABEL, "legend.fontsize": FONT_LEGEND})
    fig, (axv, axm) = plt.subplots(1, 2, figsize=(13, 5.6), dpi=450,
                                   gridspec_kw={"width_ratios": [1, 1.45], "wspace": 0.22})

    # ── Panel A: per-species violin + box of log2(B/A) ──────────────────────
    data = [w.loc[w.species == s, "log2BA"].values for s in SPECIES_ORDER]
    vp = axv.violinplot(data, showextrema=False, widths=0.85)
    for body, s in zip(vp["bodies"], SPECIES_ORDER):
        body.set_facecolor(SPECIES_COLOR[s]); body.set_alpha(0.45); body.set_edgecolor("none")
    bp = axv.boxplot(data, widths=0.18, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", lw=1.6), whiskerprops=dict(color="0.3"),
                     capprops=dict(color="0.3"), boxprops=dict(edgecolor="0.3"))
    for patch, s in zip(bp["boxes"], SPECIES_ORDER):
        patch.set_facecolor(SPECIES_COLOR[s]); patch.set_alpha(0.9)
    for i, s in enumerate(SPECIES_ORDER, start=1):
        axv.hlines(exp_log2[s], i - 0.45, i + 0.45, color=SPECIES_COLOR[s], ls="--", lw=2)
        med = float(np.median(w.loc[w.species == s, "log2BA"]))
        axv.annotate(f"obs {med:+.2f}\nexp {exp_log2[s]:+.2f}", (i, exp_log2[s]),
                     xytext=(i + 0.46, exp_log2[s]), fontsize=FONT_ANNOTATION - 3,
                     va="center", ha="left", color="0.2")
    axv.axhline(0, color="0.6", lw=0.8, zorder=0)
    axv.set_xticks(range(1, 4)); axv.set_xticklabels([SPECIES_LABEL[s] for s in SPECIES_ORDER])
    axv.set_ylabel("log2(B / A)   (human-anchored)", fontsize=FONT_AXIS_LABEL)
    axv.set_title("Ratio recovery by species", fontsize=FONT_TITLE)
    axv.set_ylim(-3, 3.2)

    # ── Panel B: LFQbench MA scatter ────────────────────────────────────────
    for s in SPECIES_ORDER:
        d = w[w.species == s]
        axm.scatter(d["log2int"], d["log2BA"], s=6, alpha=0.25,
                    color=SPECIES_COLOR[s], edgecolors="none",
                    label=f"{SPECIES_LABEL[s]}  (n={len(d):,})")
        axm.axhline(exp_log2[s], color=SPECIES_COLOR[s], ls="--", lw=1.6, alpha=0.9)
    axm.axhline(0, color="0.6", lw=0.8, zorder=0)
    axm.set_xlabel("log2 mean precursor abundance  ½·log2(A·B)", fontsize=FONT_AXIS_LABEL)
    axm.set_ylabel("log2(B / A)   (human-anchored)", fontsize=FONT_AXIS_LABEL)
    axm.set_title("LFQbench view (dashed = expected ratio)", fontsize=FONT_TITLE)
    axm.set_ylim(-3, 3.2)
    leg = axm.legend(loc="upper right", framealpha=0.9, markerscale=3, handletextpad=0.4)
    for lh in leg.legend_handles:
        lh.set_alpha(1)

    n = len(w)
    fig.suptitle(f"{args.title}\n{n:,} precursors quantified in both samples · "
                 f"expected B/A: human 1.0, yeast 3.0, E. coli 0.5",
                 fontsize=FONT_TITLE, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    print(f"wrote {args.out}  ({n:,} precursors; human anchor shift {anchor:+.3f} log2)")
    for s in SPECIES_ORDER:
        d = w[w.species == s]
        print(f"  {s:6s} n={len(d):6,d}  median log2(B/A)={np.median(d['log2BA']):+.3f}  "
              f"(expected {exp_log2[s]:+.3f})")


if __name__ == "__main__":
    main()
