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

from .ratios import (SPECIES_ORDER, SPECIES_LABEL, SPECIES_COLOR,
                     RatioConvention, A_OVER_B, EXPECTED_LOG2_BA)


def _feature(df: pd.DataFrame, level: str) -> pd.Series:
    if level in ("ion", "precursor"):
        return df["sequence_modified"].astype(str) + "@" + df["charge"].astype(str)
    if level == "modified_peptide":
        return df["sequence_modified"].astype(str)
    if level == "peptide":
        return df["sequence"].astype(str)
    raise ValueError(f"level must be ion|peptide|modified_peptide, got {level!r}")


def build_true_abundance(truth: pd.DataFrame, obs: pd.DataFrame, level: str = "ion",
                         conv: RatioConvention = A_OVER_B, n_bins: int = 10,
                         q_value_max: float = 0.01) -> pd.DataFrame:
    """Per (species, true-abundance decile): sensitivity + median ratio bias."""
    t = truth[truth["transmitted"] & (truth["truth_scope"] != "background_unknown")
              & truth["species"].isin(SPECIES_ORDER)].copy()
    t["feature"] = _feature(t, level)
    sp = t.groupby("feature")["species"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None)

    # true abundance per feature: need it in BOTH samples for the MA midpoint
    ti = t.groupby(["feature", "sample"])["truth_intensity"].sum().unstack("sample")
    ti = ti.dropna(subset=["A", "B"])
    ti = ti[(ti["A"] > 0) & (ti["B"] > 0)]          # both samples positive for the MA midpoint
    feat = pd.DataFrame(index=ti.index)
    feat["true_log2_int"] = 0.5 * np.log2(ti["A"] * ti["B"])
    feat["species"] = feat.index.map(sp)
    feat = feat.dropna(subset=["species"])

    # observations: q-filtered, per feature per sample
    o = obs.copy()
    if q_value_max is not None and "q_value" in o.columns:
        o = o[o["q_value"].fillna(1.0) <= q_value_max]
    o["feature"] = _feature(o, level)
    detected_any = set(o["feature"])
    ow = o.groupby(["feature", "sample"])["observed_intensity"].sum().unstack("sample")
    ow_both = ow.dropna(subset=["A", "B"]).copy()
    ow_both = ow_both[(ow_both["A"] > 0) & (ow_both["B"] > 0)]
    ow_both["log2_ba"] = np.log2(ow_both["B"] / ow_both["A"])
    # human-anchor the observed ratio (consistent with the ratio panels)
    h = ow_both["log2_ba"].reindex([f for f in feat.index[feat["species"] == "HUMAN"]
                                    if f in ow_both.index])
    anchor = float(h.median()) if len(h) else 0.0

    feat["detected"] = feat.index.isin(detected_any)
    feat["ratio"] = (ow_both["log2_ba"] - anchor).reindex(feat.index)

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
    return pd.DataFrame(rows)


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
