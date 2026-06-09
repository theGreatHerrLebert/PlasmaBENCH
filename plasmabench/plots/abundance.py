"""P3 — sensitivity & ratio bias vs *true* (blueprint) abundance.

The cross-version-correct abundance analysis: bin precursors by their **ground-truth**
intensity from the blueprint (identical for every engine/version, so the x-axis is
shared and comparable — unlike each engine's own reported abundance), then per
species per bin report:

  * **detection sensitivity** = fraction of transmitted truth precursors the engine
    identified (the real LOD/sensitivity curve), and
  * **ratio bias** = median observed log2(A/B) among precursors quantified in both
    samples, vs the expected ratio.

True abundance of a precursor = ``½·log2(truth_A · truth_B)`` (the MA midpoint, from
the blueprint). Supports ion / peptide / modified_peptide (keys that match cleanly
between truth and observations).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .ratios import (SPECIES_ORDER, SPECIES_LABEL, SPECIES_COLOR, Q_COLUMNS,
                     RatioConvention, A_OVER_B, EXPECTED_LOG2_BA, build_ratio_table)

ABUND_COLUMNS = ["species", "x", "n_truth", "n_detected", "sensitivity",
                 "median_ratio", "q25", "q75"]


def _feature(df: pd.DataFrame, level: str) -> pd.Series:
    """Feature key matching ratios.build_ratio_table exactly (so ratios join)."""
    if level in ("ion", "precursor"):
        return df["sequence_modified"].astype(str) + df["charge"].astype(str)
    if level == "modified_peptide":
        return df["sequence_modified"].astype(str)
    if level == "peptide":
        return df["sequence"].astype(str)
    raise ValueError(f"level must be ion|peptide|modified_peptide, got {level!r}")


def _found_features(obs: pd.DataFrame, level: str, q_value_max: float) -> set:
    """Features identified in ≥1 sample, using the SAME 4-level q filter as
    build_ratio_table (null q-values permitted)."""
    o = obs
    if q_value_max is not None:
        for qc in Q_COLUMNS:
            if qc in o.columns:
                o = o[o[qc].isna() | (o[qc] <= q_value_max)]
    o = o[o["observed_intensity"] > 0]
    return set(_feature(o, level))


def build_true_abundance(truth: pd.DataFrame, obs: pd.DataFrame, level: str = "ion",
                         conv: RatioConvention = A_OVER_B, n_bins: int = 10,
                         q_value_max: float = 0.01) -> pd.DataFrame:
    """Per (species, true-abundance decile): detection sensitivity + median ratio bias.

    Sensitivity = P(identified in ≥1 sample | transmitted in ≥1 sample), over the
    findable truth universe. The bias ratio is taken from build_ratio_table, so it
    is identical to the P1 ratio panels (same filters, aggregation, human anchor).
    """
    for df_, nm in ((truth, "truth"), (obs, "observations")):
        if df_["experiment"].nunique() > 1:
            raise ValueError(f"{nm} spans >1 experiment; call build_true_abundance "
                             "per experiment (feature keys are not experiment-scoped)")

    # observed ratio per feature — reuse the P1 ratio table for exact consistency
    wide, _ = build_ratio_table(obs, level=level, q_value_max=q_value_max, human_anchor=True)
    ratio_by_feat = wide.set_index("feature")["log2_ba"]
    found = _found_features(obs, level, q_value_max)

    # truth universe: positive intensity in BOTH samples (for the MA midpoint),
    # findable = transmitted in ≥1 sample (NOT both — that would drop weak features)
    t = truth[(truth["truth_scope"] != "background_unknown")
              & truth["species"].isin(SPECIES_ORDER) & (truth["truth_intensity"] > 0)].copy()
    t["feature"] = _feature(t, level)
    sp = t.groupby("feature")["species"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None)
    transmitted_any = t.groupby("feature")["transmitted"].any()
    ti = t.groupby(["feature", "sample"])["truth_intensity"].sum().unstack("sample")
    ti = ti.dropna(subset=["A", "B"])
    ti = ti[(ti["A"] > 0) & (ti["B"] > 0)]

    feat = pd.DataFrame(index=ti.index)
    feat["true_log2_int"] = 0.5 * np.log2(ti["A"] * ti["B"])
    feat["species"] = feat.index.map(sp)
    feat["eligible"] = feat.index.map(transmitted_any).fillna(False)
    feat = feat.dropna(subset=["species"])
    feat = feat[feat["eligible"]]                 # findable universe = the sensitivity denominator
    feat["detected"] = feat.index.isin(found)
    feat["ratio"] = ratio_by_feat.reindex(feat.index)   # anchored canonical log2(B/A); NaN if not quantified in both

    rows = []
    for s in SPECIES_ORDER:
        d = feat[feat["species"] == s]
        if len(d) < n_bins:
            continue
        abin = pd.qcut(d["true_log2_int"], n_bins, duplicates="drop")
        for _, g in d.groupby(abin, observed=True):
            ratios = conv.to_display(g["ratio"].dropna())
            rows.append({
                "species": s, "x": float(g["true_log2_int"].median()),
                "n_truth": int(len(g)), "n_detected": int(g["detected"].sum()),
                "sensitivity": float(g["detected"].mean()),
                "median_ratio": float(ratios.median()) if len(ratios) else np.nan,
                "q25": float(ratios.quantile(0.25)) if len(ratios) else np.nan,
                "q75": float(ratios.quantile(0.75)) if len(ratios) else np.nan,
            })
    return pd.DataFrame(rows, columns=ABUND_COLUMNS)


def plot_true_abundance(tables: dict, out: Path, conv: RatioConvention = A_OVER_B,
                        level: str = "ion", y_halfspan: float = 1.2,
                        title: str = "Sensitivity & ratio bias vs TRUE abundance") -> pd.DataFrame:
    """2 rows (sensitivity, ratio bias) × 3 species; one line per label in ``tables``."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"axes.labelsize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10})
    palette = ["#B07AA1", "#4E79A7", "#59A14F", "#E15759"]
    fig, axes = plt.subplots(2, len(SPECIES_ORDER), figsize=(5.0 * len(SPECIES_ORDER), 8.4),
                             dpi=300, sharex="col")

    for j, s in enumerate(SPECIES_ORDER):
        ax_sens, ax_bias = axes[0][j], axes[1][j]
        exp = conv.to_display(EXPECTED_LOG2_BA[s])
        for i, (lbl, tb) in enumerate(tables.items()):
            d = tb[tb["species"] == s].sort_values("x")
            if d.empty:
                continue
            c = palette[i % len(palette)]
            ax_sens.plot(d["x"], 100 * d["sensitivity"], "-o", color=c, ms=4, lw=1.8, label=lbl)
            db = d.dropna(subset=["median_ratio"])   # only bins with quantified-in-both precursors
            ax_bias.fill_between(db["x"], db["q25"], db["q75"], color=c, alpha=0.12, lw=0)
            ax_bias.plot(db["x"], db["median_ratio"], "-o", color=c, ms=4, lw=1.8, label=lbl)
        ax_sens.set_title(SPECIES_LABEL[s]); ax_sens.set_ylim(0, 102)
        ax_sens.grid(axis="y", ls=":", alpha=0.4); ax_sens.legend(fontsize=9, loc="lower right")
        ax_bias.axhline(exp, color="0.4", ls="--", lw=1.4, label=f"expected {exp:+.2f}")
        ax_bias.set_ylim(exp - y_halfspan, exp + y_halfspan)
        ax_bias.grid(axis="y", ls=":", alpha=0.4); ax_bias.legend(fontsize=9, loc="best")
        ax_bias.set_xlabel("log2 TRUE mean abundance  ½·log2(tA·tB)")
    axes[0][0].set_ylabel("detection sensitivity %")
    axes[1][0].set_ylabel(conv.ylabel + " (median per decile)")
    fig.suptitle(f"{title} — {level} level", fontsize=14, y=1.0)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return pd.concat([t.assign(label=lbl) for lbl, t in tables.items()], ignore_index=True)


__all__ = ["build_true_abundance", "plot_true_abundance"]
