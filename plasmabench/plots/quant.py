"""P6 — ground-truth quant correlation (precursor level).

The panel the LFQ-Bench paper structurally cannot produce: simulated truth
intensity vs engine-recovered quantity, per species, against exact truth.

Metric (decided in plan.md, do not silently change):
  * **primary:** Spearman(log10 truth, log10 observed) per ``(sample, species)`` —
    rank-based, so it is immune to the absolute-unit mismatch between blueprint
    ``events × relative_abundance`` and the engine's reported quantity.
  * **secondary:** log-Pearson after a per-run scale normalization (centre each
    run's log-observed on its median) — reported, not headline.

Joins ``truth`` (transmitted only — the findable universe) to ``observations``
(q-value filtered) on ``(experiment, sample, run_id, sequence_modified, charge)``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# timsim-bench house style, kept consistent with scripts/plot_stage1_ratios.py
SPECIES_COLOR = {"HUMAN": "#DD8452", "YEAST": "#4C72B0", "ECOLI": "#55A868"}  # plasma orange, yeast blue
SPECIES_ORDER = ["ECOLI", "HUMAN", "YEAST"]
SPECIES_LABEL = {"HUMAN": "Human (plasma)", "YEAST": "Yeast", "ECOLI": "E. coli"}

JOIN_KEYS = ["experiment", "sample", "run_id", "sequence_modified", "charge"]


def build_quant_table(truth: pd.DataFrame, obs: pd.DataFrame,
                      q_value_max: float = 0.01,
                      transmitted_only: bool = True) -> pd.DataFrame:
    """Inner-join transmitted truth to q-filtered observations → per-precursor points.

    Species is taken from **truth** (the oracle); rows where truth species is
    ambiguous (None) are dropped — they are unscorable, never counted as errors.
    """
    t = truth.copy()
    if transmitted_only:
        t = t[t["transmitted"]]
    t = t[t["species"].isin(SPECIES_ORDER)]
    t = t[t["truth_intensity"] > 0]

    o = obs.copy()
    if q_value_max is not None and "q_value" in o.columns:
        o = o[o["q_value"].fillna(1.0) <= q_value_max]
    o = o[o["observed_intensity"] > 0]

    merged = t.merge(
        o[JOIN_KEYS + ["observed_intensity"]],
        on=JOIN_KEYS, how="inner", validate="one_to_one",
    )
    merged["log_truth"] = np.log10(merged["truth_intensity"])
    merged["log_obs"] = np.log10(merged["observed_intensity"])
    return merged


def correlations(merged: pd.DataFrame) -> pd.DataFrame:
    """Per ``(sample, species)`` Spearman (primary) + scale-normalized log-Pearson."""
    rows = []
    for (sample, species), g in merged.groupby(["sample", "species"], sort=False):
        if len(g) < 3:
            continue
        rho = spearmanr(g["log_truth"], g["log_obs"]).statistic
        # per-run scale normalization, then log-Pearson
        gg = g.copy()
        gg["log_obs_norm"] = gg["log_obs"] - gg.groupby("run_id")["log_obs"].transform("median")
        gg["log_truth_norm"] = gg["log_truth"] - gg.groupby("run_id")["log_truth"].transform("median")
        r = pearsonr(gg["log_truth_norm"], gg["log_obs_norm"]).statistic
        rows.append({"sample": sample, "species": species, "n": len(g),
                     "spearman": rho, "log_pearson_scaled": r})
    return pd.DataFrame(rows).sort_values(["sample", "species"]).reset_index(drop=True)


def plot_quant_correlation(merged: pd.DataFrame, corrs: pd.DataFrame, out: Path,
                           title: str = "PlasmaBENCH P6 — quant correlation (precursor)") -> None:
    """One scatter per sample: log10 truth vs log10 observed, coloured by species."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11})
    samples = list(dict.fromkeys(merged["sample"]))
    fig, axes = plt.subplots(1, len(samples), figsize=(6.2 * len(samples), 5.6), dpi=300, squeeze=False)
    corr_lookup = {(r["sample"], r["species"]): r for _, r in corrs.iterrows()}

    for ax, sample in zip(axes[0], samples):
        d = merged[merged["sample"] == sample]
        for sp in SPECIES_ORDER:
            ds = d[d["species"] == sp]
            if ds.empty:
                continue
            c = corr_lookup.get((sample, sp))
            lbl = f"{SPECIES_LABEL[sp]} (n={len(ds):,}"
            lbl += f", ρ={c['spearman']:.2f})" if c is not None else ")"
            ax.scatter(ds["log_truth"], ds["log_obs"], s=6, alpha=0.25,
                       color=SPECIES_COLOR[sp], edgecolors="none", label=lbl)
        lo = float(np.nanmin([d["log_truth"].min(), d["log_obs"].min()]))
        hi = float(np.nanmax([d["log_truth"].max(), d["log_obs"].max()]))
        ax.plot([lo, hi], [lo, hi], color="0.5", ls="--", lw=1, zorder=0)
        ax.set_xlabel("log10 simulated truth intensity")
        ax.set_ylabel("log10 observed quantity")
        ax.set_title(f"Sample {sample}")
        leg = ax.legend(loc="upper left", framealpha=0.9, markerscale=3, handletextpad=0.4)
        for lh in leg.legend_handles:
            lh.set_alpha(1)

    fig.suptitle(title, fontsize=14, y=1.0)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


__all__ = ["build_quant_table", "correlations", "plot_quant_correlation", "JOIN_KEYS"]
