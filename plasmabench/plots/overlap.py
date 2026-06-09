"""Cross-engine ID overlap — 2-way Venn at precursor / peptide / protein level.

How much do two DIA software / versions agree on *what* they identify? Reads each
engine's report, builds the identified-ID set at each level (q-filtered, pooled
across the runs in the report), and draws a 2-circle Venn per level with exact
region counts. This is the cross-engine *difference* view the benchmark cares about
(see CLAUDE.md "cross-engine integrity") and mirrors the paper's overlap panels.

Circles are schematic (fixed radius); the annotated counts are exact.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .loader import read_diann_minimal

LEVELS = ["precursor", "peptide", "protein"]
LEVEL_LABEL = {"precursor": "Precursor (seq+charge)", "peptide": "Peptide (stripped)",
               "protein": "Protein (accession)"}


def id_sets(report_path: str | Path, q_value_max: float = 0.01) -> dict[str, set]:
    """Identified-ID sets at each level for one report, pooled across its runs."""
    df = read_diann_minimal(report_path, q_value_max)
    prec = set(zip(df["sequence_modified"], df["charge"]))
    pep = set(df["sequence"].dropna())
    prot: set[str] = set()
    for g in df["protein_group"].dropna():
        for tok in str(g).split(";"):
            tok = tok.strip()
            if tok:
                prot.add(tok)
    return {"precursor": prec, "peptide": pep, "protein": prot}


def overlap_counts(a: set, b: set) -> dict:
    """``{only_a, shared, only_b, jaccard}`` for two sets."""
    shared = len(a & b)
    union = len(a | b)
    return {"only_a": len(a - b), "shared": shared, "only_b": len(b - a),
            "jaccard": (shared / union) if union else 0.0}


def plot_overlap_venn(sets_a: dict, sets_b: dict, label_a: str, label_b: str, out: Path,
                      color_a: str = "#B07AA1", color_b: str = "#4E79A7",
                      title: str = "ID overlap") -> "pd.DataFrame":
    """3-panel (precursor/peptide/protein) 2-circle Venn; returns the counts table."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pandas as pd

    mpl.rcParams.update({"font.size": 12})
    fig, axes = plt.subplots(1, len(LEVELS), figsize=(5.2 * len(LEVELS), 5.0), dpi=300)
    rows = []
    for ax, level in zip(axes, LEVELS):
        a, b = sets_a[level], sets_b[level]
        c = overlap_counts(a, b)
        rows.append({"level": level, "n_" + label_a: len(a), "n_" + label_b: len(b),
                     "only_" + label_a: c["only_a"], "shared": c["shared"],
                     "only_" + label_b: c["only_b"], "jaccard": round(c["jaccard"], 4)})
        r = 0.62
        ax.add_patch(plt.Circle((-0.30, 0), r, color=color_a, alpha=0.45, ec="none"))
        ax.add_patch(plt.Circle((0.30, 0), r, color=color_b, alpha=0.45, ec="none"))
        ax.text(-0.62, 0, f"{c['only_a']:,}", ha="center", va="center", fontsize=14)
        ax.text(0.0, 0, f"{c['shared']:,}", ha="center", va="center", fontsize=15, fontweight="bold")
        ax.text(0.62, 0, f"{c['only_b']:,}", ha="center", va="center", fontsize=14)
        ax.text(-0.45, r + 0.07, label_a, ha="center", va="bottom", color=color_a, fontweight="bold")
        ax.text(0.45, r + 0.07, label_b, ha="center", va="bottom", color=color_b, fontweight="bold")
        ax.set_title(f"{LEVEL_LABEL[level]}\n{len(a):,} vs {len(b):,}  ·  "
                     f"{100 * c['jaccard']:.1f}% shared", fontsize=12)
        ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.0, 1.2)
        ax.set_aspect("equal"); ax.axis("off")

    fig.suptitle(title, fontsize=15, y=1.02)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(rows)


__all__ = ["id_sets", "overlap_counts", "plot_overlap_venn", "LEVELS"]
