"""P1 / P4 — three-species A/B ratio panels (violin + LFQbench MA scatter).

The ported core of the PYE paper's Fig_5A / Fig_5B_C / Fig_6: per-species
distributions of the Sample-A/Sample-B fold change against the expected PYE1
ratios. Engine-internal (no blueprint needed) — it reads the ``observations``
table only.

Decisions baked in (see plan.md):
  * **Protein-level (PG.MaxLFQ) is primary**; ``level`` also offers the lower tiers
    of the ID hierarchy on ``Precursor.Normalised``:
      - ``protein``          — `protein_group`, PG.MaxLFQ (dedup per run)
      - ``peptide``          — stripped sequence; sum precursors over charges + modforms
      - ``modified_peptide`` — peptidoform `sequence_modified`; sum over charges
      - ``ion`` (=``precursor``) — `(sequence_modified, charge)`; the precursor itself
  * **One display API.** The canonical quantity is ``log2(B/A)`` (matches
    ``EXPECTED_RATIO_B_OVER_A``); a :class:`RatioConvention` flips it for display.
    Data *and* expected reference lines go through the **same** ``to_display`` —
    they cannot drift. Default convention is ``A/B`` (the paper's axis: human 0,
    yeast −log2 3, E. coli +1).
  * **Aggregation across replicates is a mode** {linear_mean, geometric_mean,
    median}; linear_mean reproduces the paper, and is a no-op for 1 run/sample.
  * Paper-parity filters: 4-level DIA-NN FDR (≤ q), ``_CONTA`` removal,
    ≥N non-NA replicates per condition, ambiguous-species excluded.
  * **Human-anchoring** (shift human median → expected) is our robustness
    deviation, on by default; toggle with ``human_anchor``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from plasmabench.species import EXPECTED_RATIO_B_OVER_A

SPECIES_COLOR = {"HUMAN": "#DD8452", "YEAST": "#4C72B0", "ECOLI": "#55A868"}  # plasma orange, yeast blue
SPECIES_ORDER = ["ECOLI", "HUMAN", "YEAST"]
SPECIES_LABEL = {"HUMAN": "Human (plasma)", "YEAST": "Yeast", "ECOLI": "E. coli"}
Q_COLUMNS = ["q_value", "lib_q_value", "pg_q_value", "lib_pg_q_value"]


@dataclass(frozen=True)
class RatioConvention:
    """How to display the canonical ``log2(B/A)``.  ``display = sign · canonical``."""
    name: str
    sign: int
    ylabel: str

    def to_display(self, canonical):
        return self.sign * canonical


A_OVER_B = RatioConvention("A/B", -1, "log2(A / B)")   # paper: human 0, yeast −log2 3, ecoli +1
B_OVER_A = RatioConvention("B/A", +1, "log2(B / A)")   # oracle-native: human 0, yeast +log2 3, ecoli −1

EXPECTED_LOG2_BA = {s: float(np.log2(EXPECTED_RATIO_B_OVER_A[s])) for s in SPECIES_ORDER}


def _agg_func(aggregation: str):
    if aggregation == "linear_mean":
        return lambda x: np.mean(x)
    if aggregation == "median":
        return lambda x: np.median(x)
    if aggregation == "geometric_mean":
        return lambda x: float(np.exp(np.mean(np.log(x))))  # x already > 0
    raise ValueError(f"aggregation must be linear_mean|geometric_mean|median, got {aggregation!r}")


def build_ratio_table(obs: pd.DataFrame, level: str = "protein",
                      aggregation: str = "linear_mean", min_replicates: int = 1,
                      q_value_max: float = 0.01, human_anchor: bool = True,
                      protein_quant: str = "maxlfq"):
    """Per-feature ``log2(B/A)`` (canonical) + mean log2 abundance, per species.

    Returns ``(wide, anchor)`` where ``wide`` has one row per feature with
    ``species, log2_ba, log2_int`` and ``anchor`` is the human-median shift applied
    (0.0 if ``human_anchor`` is False).

    ``protein_quant`` (level=='protein' only): ``"maxlfq"`` uses DIA-NN's PG.MaxLFQ
    (carries the engine's cross-run normalization); ``"raw_sum"`` sums the raw
    ``observed_intensity`` over the protein group — the normalization-free protein
    quantity used when the benchmark is run on the raw quantity axis (so 1.8 and 2.5
    agree at protein level too). See docs/finding-diann-1.8-vs-2.5-quant.md.
    """
    if protein_quant not in ("maxlfq", "raw_sum"):
        raise ValueError(f"protein_quant must be maxlfq|raw_sum, got {protein_quant!r}")
    df = obs.copy()

    # 4-level FDR filter (apply each q-column that is present and non-null).
    for qc in Q_COLUMNS:
        if qc in df.columns:
            df = df[df[qc].isna() | (df[qc] <= q_value_max)]
    # contaminants & ambiguous species out (assign_species already drops _CONTA).
    if "protein_names" in df.columns:
        df = df[~df["protein_names"].fillna("").str.contains("_CONTA")]
    df = df[df["species"].isin(SPECIES_ORDER)]

    # Feature key per level (the unit whose A/B ratio we plot).
    if level == "protein":
        if protein_quant == "maxlfq":
            df = df.dropna(subset=["protein_group", "pg_maxlfq"])
            df = df[df["pg_maxlfq"] > 0]
        else:  # raw_sum: aggregate raw precursor observed_intensity over the group.
            # Purpose is a normalization-free protein quantity for raw-vs-raw cross-engine
            # comparison — NOT to reproduce MaxLFQ. Like maxlfq it sums only q-passing
            # precursors, so both modes share the same precursor inclusion; they differ
            # only in aggregation (raw sum vs MaxLFQ's normalized group quantity).
            df = df.dropna(subset=["protein_group", "observed_intensity"])
            df = df[df["observed_intensity"] > 0]
        df["feature"] = df["protein_group"]
    elif level in ("ion", "precursor"):  # precursor ion = modified sequence + charge
        df = df.dropna(subset=["observed_intensity"])
        df = df[df["observed_intensity"] > 0]
        df["feature"] = df["sequence_modified"] + df["charge"].astype(str)
    elif level == "modified_peptide":    # peptidoform: collapse charges
        df = df.dropna(subset=["observed_intensity"])
        df = df[df["observed_intensity"] > 0]
        df["feature"] = df["sequence_modified"]
    elif level == "peptide":             # stripped sequence: collapse charges + modforms
        df = df.dropna(subset=["observed_intensity", "sequence"])
        df = df[df["observed_intensity"] > 0]
        df["feature"] = df["sequence"]
    else:
        raise ValueError(f"level must be protein|peptide|modified_peptide|ion, got {level!r}")

    # feature → species (drop features whose species is ambiguous: cross-species leakage guard)
    feat_species = df.groupby("feature")["species"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna()

    # Within a (run, feature): protein uses PG.MaxLFQ (one per group) or, in raw_sum
    # mode, the sum of raw precursor intensities; ion is already one row per precursor;
    # peptide/peptidoform SUM precursor intensities.
    if level == "protein" and protein_quant == "maxlfq":
        per_run = df.drop_duplicates(["run_id", "feature"])[
            ["sample", "run_id", "feature", "pg_maxlfq"]].rename(columns={"pg_maxlfq": "quantity"})
    elif level == "protein":  # raw_sum
        per_run = df.groupby(["sample", "run_id", "feature"], as_index=False)["observed_intensity"] \
            .sum().rename(columns={"observed_intensity": "quantity"})
    elif level in ("ion", "precursor"):
        per_run = df[["sample", "run_id", "feature", "observed_intensity"]] \
            .rename(columns={"observed_intensity": "quantity"})
    else:  # peptide / modified_peptide → sum precursor intensities within the run
        per_run = df.groupby(["sample", "run_id", "feature"], as_index=False)["observed_intensity"] \
            .sum().rename(columns={"observed_intensity": "quantity"})

    # aggregate replicates within each (feature, sample)
    af = _agg_func(aggregation)
    agg = per_run.groupby(["feature", "sample"])["quantity"].agg(qty=af, count="count").reset_index()

    wide = agg.pivot(index="feature", columns="sample", values=["qty", "count"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.dropna(subset=["qty_A", "qty_B"])
    wide = wide[(wide["count_A"] >= min_replicates) & (wide["count_B"] >= min_replicates)]
    wide["species"] = wide.index.map(feat_species)
    wide = wide.dropna(subset=["species"])

    wide["log2_ba"] = np.log2(wide["qty_B"] / wide["qty_A"])
    wide["log2_int"] = 0.5 * np.log2(wide["qty_A"] * wide["qty_B"])

    anchor = 0.0
    if human_anchor:
        h = wide.loc[wide["species"] == "HUMAN", "log2_ba"]
        if len(h):
            anchor = float(h.median())
            wide["log2_ba"] = wide["log2_ba"] - anchor
    return wide.reset_index(), anchor


def ratio_summary(wide: pd.DataFrame, conv: RatioConvention) -> pd.DataFrame:
    """Per-species observed vs expected median ratio, in the display convention."""
    rows = []
    for s in SPECIES_ORDER:
        d = wide[wide["species"] == s]
        if d.empty:
            continue
        rows.append({
            "species": s, "n": len(d),
            "observed_log2": conv.to_display(float(d["log2_ba"].median())),
            "expected_log2": conv.to_display(EXPECTED_LOG2_BA[s]),
        })
    return pd.DataFrame(rows)


def plot_ratio_panels(wide: pd.DataFrame, out: Path, conv: RatioConvention = A_OVER_B,
                      level: str = "protein", anchor: float = 0.0,
                      title: str = "PlasmaBENCH P1 — three-species ratio") -> None:
    """Violin (per species) + LFQbench MA scatter, in the chosen display convention."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12})
    fig, (axv, axm) = plt.subplots(1, 2, figsize=(13, 5.6), dpi=300,
                                   gridspec_kw={"width_ratios": [1, 1.45], "wspace": 0.22})
    disp = {s: conv.to_display(wide.loc[wide.species == s, "log2_ba"].values) for s in SPECIES_ORDER}
    exp_disp = {s: conv.to_display(EXPECTED_LOG2_BA[s]) for s in SPECIES_ORDER}

    data = [disp[s] for s in SPECIES_ORDER]
    vp = axv.violinplot(data, showextrema=False, widths=0.85)
    for body, s in zip(vp["bodies"], SPECIES_ORDER):
        body.set_facecolor(SPECIES_COLOR[s]); body.set_alpha(0.45); body.set_edgecolor("none")
    bp = axv.boxplot(data, widths=0.18, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", lw=1.6), whiskerprops=dict(color="0.3"),
                     capprops=dict(color="0.3"), boxprops=dict(edgecolor="0.3"))
    for patch, s in zip(bp["boxes"], SPECIES_ORDER):
        patch.set_facecolor(SPECIES_COLOR[s]); patch.set_alpha(0.9)
    for i, s in enumerate(SPECIES_ORDER, start=1):
        axv.hlines(exp_disp[s], i - 0.45, i + 0.45, color=SPECIES_COLOR[s], ls="--", lw=2)
        med = float(np.median(disp[s]))
        axv.annotate(f"obs {med:+.2f}\nexp {exp_disp[s]:+.2f}", (i, exp_disp[s]),
                     xytext=(i + 0.46, exp_disp[s]), fontsize=11, va="center", ha="left", color="0.2")
    axv.axhline(0, color="0.6", lw=0.8, zorder=0)
    axv.set_xticks(range(1, 4)); axv.set_xticklabels([SPECIES_LABEL[s] for s in SPECIES_ORDER])
    axv.set_ylabel(conv.ylabel + ("   (human-anchored)" if anchor else ""))
    axv.set_title(f"Ratio recovery by species ({level})")
    axv.set_ylim(-3, 3.2)

    # Protein has far fewer points than peptide/ion — give it larger, more
    # opaque markers so the spread is legible; keep the dense levels small.
    pt_size, pt_alpha = (28, 0.42) if level == "protein" else (6, 0.25)
    for s in SPECIES_ORDER:
        d = wide[wide.species == s]
        axm.scatter(d["log2_int"], conv.to_display(d["log2_ba"]), s=pt_size, alpha=pt_alpha,
                    color=SPECIES_COLOR[s], edgecolors="none",
                    label=f"{SPECIES_LABEL[s]}  (n={len(d):,})")
        axm.axhline(exp_disp[s], color=SPECIES_COLOR[s], ls="--", lw=1.6, alpha=0.9)
    axm.axhline(0, color="0.6", lw=0.8, zorder=0)
    axm.set_xlabel("log2 mean abundance  ½·log2(A·B)")
    axm.set_ylabel(conv.ylabel)
    axm.set_title("LFQbench view (dashed = expected ratio)")
    axm.set_ylim(-3, 3.2)
    leg = axm.legend(loc="upper right", framealpha=0.9,
                     markerscale=2 if level == "protein" else 3, handletextpad=0.4)
    for lh in leg.legend_handles:
        lh.set_alpha(1)

    n = len(wide)
    fig.suptitle(f"{title}\n{n:,} {level}s quantified in both samples · "
                 f"expected {conv.name}: " + ", ".join(
                     f"{SPECIES_LABEL[s].split()[0].lower()} {2**exp_disp[s]:.2g}" for s in SPECIES_ORDER),
                 fontsize=14, y=0.99)
    fig.subplots_adjust(top=0.86, left=0.07, right=0.98, bottom=0.11)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def bias_by_abundance(wide: pd.DataFrame, conv: RatioConvention = A_OVER_B,
                      n_bins: int = 10) -> pd.DataFrame:
    """Per-species, per-abundance-decile median ratio (display convention) + IQR.

    Bins are equal-count quantiles of ``log2_int`` *within* each species (their
    abundance ranges differ). Reveals abundance-dependent ratio bias.
    """
    rows = []
    for sp in SPECIES_ORDER:
        d = wide[wide["species"] == sp]
        if len(d) < n_bins:
            continue
        abin = pd.qcut(d["log2_int"], n_bins, duplicates="drop")
        for _, g in d.groupby(abin, observed=True):
            disp = conv.to_display(g["log2_ba"])
            rows.append({"species": sp, "x": float(g["log2_int"].median()),
                         "median": float(disp.median()),
                         "q25": float(disp.quantile(0.25)), "q75": float(disp.quantile(0.75)),
                         "n": int(len(g))})
    return pd.DataFrame(rows, columns=["species", "x", "median", "q25", "q75", "n"])


def plot_ratio_bias(wides: dict, out: Path, conv: RatioConvention = A_OVER_B,
                    level: str = "protein", n_bins: int = 10, y_halfspan: float = 1.0,
                    title: str = "Ratio bias vs abundance") -> "pd.DataFrame":
    """One panel per species: median ratio vs abundance decile, overlaying each
    label in ``wides`` (e.g. {"1.8": wide18, "2.5": wide25}). Dashed = expected.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11})
    palette = ["#B07AA1", "#4E79A7", "#59A14F", "#E15759"]
    fig, axes = plt.subplots(1, len(SPECIES_ORDER), figsize=(5.2 * len(SPECIES_ORDER), 5.0), dpi=300)
    tables = {lbl: bias_by_abundance(w, conv, n_bins) for lbl, w in wides.items()}

    for ax, sp in zip(axes, SPECIES_ORDER):
        exp = conv.to_display(EXPECTED_LOG2_BA[sp])
        ax.axhline(exp, color="0.4", ls="--", lw=1.4, label=f"expected {exp:+.2f}")
        for i, (lbl, tb) in enumerate(tables.items()):
            d = tb[tb["species"] == sp]
            if d.empty:
                continue
            c = palette[i % len(palette)]
            ax.fill_between(d["x"], d["q25"], d["q75"], color=c, alpha=0.12, lw=0)
            ax.plot(d["x"], d["median"], "-o", color=c, ms=4, lw=1.8, label=lbl)
        ax.set_title(f"{SPECIES_LABEL[sp]}")
        ax.set_xlabel("log2 mean abundance  ½·log2(A·B)")
        ax.set_ylim(exp - y_halfspan, exp + y_halfspan)
        ax.grid(axis="y", ls=":", alpha=0.4)
        ax.legend(fontsize=10, loc="best")
    axes[0].set_ylabel(conv.ylabel + " (median per decile)")
    fig.suptitle(f"{title} — {level} level", fontsize=14, y=1.02)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return pd.concat([t.assign(label=lbl) for lbl, t in tables.items()], ignore_index=True)


__all__ = [
    "RatioConvention", "A_OVER_B", "B_OVER_A", "EXPECTED_LOG2_BA",
    "build_ratio_table", "ratio_summary", "plot_ratio_panels",
    "bias_by_abundance", "plot_ratio_bias",
]
