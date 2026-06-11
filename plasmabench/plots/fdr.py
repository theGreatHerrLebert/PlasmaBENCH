"""P5 — simulated-spike-in FDR / TPR (the panel the paper structurally cannot do).

Compares the engine's *found set* to the TimSim *truth set* from the blueprint —
true positives, false positives, and recall measured against exact ground truth,
not the engine's self-reported q-value.

Adapted port of timsim-bench ``benchmark/fpr_tpr_plots.py`` (``calculate_fdr`` /
``calculate_tpr`` / the group-credit protein rule), source commit ``08cfe24``,
with the corrections required here (see plan.md):

  * **Scoped per ``(experiment, sample)``** — truth and found are matched within a
    sample, never unioned across A/B.
  * **Set comparison, not an inner join.** FP = found − truth, FN = truth − found.
  * **`transmitted` truth only** is the recall denominator (untransmitted precursors
    are unfindable and must not count as false negatives).
  * **Restricted to scorable species** (``truth_scope == 'simulated'``). In the
    YE-only-over-real-plasma condition human is ``background_unknown`` and is
    dropped from scoring — hence *simulated-spike-in* FDR, not "real" FDR. A found
    item with ambiguous species (None) is **unscorable**, never a false positive.
  * **Protein group-credit scores the reported group as one unit** — any member ∩
    truth ⇒ TP group (primary); strict (all members ∈ truth) is a reported bound.
    Non-truth candidates are never expanded into the TP protein set.

Protein identity note: our `from_findings` blueprint stores **entry names**
(``CFAB_HUMAN``), which match DIA-NN ``Protein.Names`` — so protein membership is
tested on the names field (``protein_names``), not the accession-based group id.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

SPECIES_ORDER = ["ECOLI", "HUMAN", "YEAST"]
# We score every finding EXCEPT those confidently attributable to an unscorable
# real background (human in the YE-only-over-real-plasma condition). Ambiguous /
# unrecognized findings are NOT excluded — absent from truth they are real FPs.
EXCLUDE_SCOPES = ("background_unknown",)


def calculate_fdr(found: set, truth: set) -> float:
    """FP / (FP + TP) × 100 — empirical FDR of a found set vs ground truth."""
    fp = len(found - truth)
    tp = len(found & truth)
    return 0.0 if (fp + tp) == 0 else round(100.0 * fp / (fp + tp), 3)


def calculate_tpr(found: set, truth: set) -> float:
    """TP / (TP + FN) × 100 — recall of a found set vs ground truth."""
    tp = len(found & truth)
    fn = len(truth - found)
    return 0.0 if (tp + fn) == 0 else round(100.0 * tp / (tp + fn), 3)


def _tokens(field) -> list[str]:
    """Protein tokens from a `;`/`,`-delimited names field.

    Does NOT split on `|` — a `sp|P12345|NAME_HUMAN` header is one protein, not
    three tokens. (Our blueprint/DIA-NN names are bare entry names anyway.)
    """
    if not isinstance(field, str) or not field:
        return []
    out = []
    for tok in re.split(r"[;,]\s*", field.strip()):
        tok = tok.strip().strip("{}[]()")
        if tok:
            out.append(tok)
    return out


def _level_keys(df: pd.DataFrame, level: str):
    if level in ("ion", "precursor"):
        return set(zip(df["sequence_modified"], df["charge"]))
    if level == "modified_peptide":
        return set(df["sequence_modified"])
    if level == "peptide":
        return set(df["sequence"])
    raise ValueError(f"set-level must be ion|peptide|modified_peptide, got {level!r}")


def _scorable_truth(truth: pd.DataFrame, experiment: str, sample: str, exclude) -> pd.DataFrame:
    # transmitted simulated truth, minus the unscorable real background. Ambiguous
    # species (None → 'unscorable') are kept: they are still real simulated peptides.
    return truth[(truth["experiment"] == experiment) & (truth["sample"] == sample)
                 & truth["transmitted"] & ~truth["truth_scope"].isin(exclude)]


def _scorable_obs(obs: pd.DataFrame, experiment: str, sample: str, exclude, q_value_max) -> pd.DataFrame:
    o = obs[(obs["experiment"] == experiment) & (obs["sample"] == sample)
            & ~obs["truth_scope"].isin(exclude)]
    if q_value_max is not None and "q_value" in o.columns:
        o = o[o["q_value"].fillna(1.0) <= q_value_max]
    return o


def blank_keysets(blank_report, q_value_max: float = 0.01,
                  levels=("ion", "peptide")) -> dict:
    """Level→set of IDs DIA-NN finds in the BLANK run (the noise source).

    These are contamination from the blank we sample noise from, not simulator
    output — subtracted from the FP side of the SIM's FDR so they don't inflate it.
    Accepts a DIA-NN report path or a pre-built observations DataFrame.
    """
    if isinstance(blank_report, (str, Path)):
        from .loader import read_diann_minimal
        df = read_diann_minimal(blank_report, q_value_max)
    else:
        df = blank_report.copy()
        if q_value_max is not None and "q_value" in df.columns:
            df = df[df["q_value"].fillna(1.0) <= q_value_max]
    return {lvl: _level_keys(df, lvl) for lvl in levels}


def _score_set_level(truth, obs, experiment, sample, level, q_value_max, exclude,
                     blank_keys=None) -> dict:
    t_keys = _level_keys(_scorable_truth(truth, experiment, sample, exclude), level)
    f_keys = _level_keys(_scorable_obs(obs, experiment, sample, exclude, q_value_max), level)
    # Subtract blank-explained IDs from the FP side only (never remove a TP that
    # also appears in the blank): found_adj = found − (blank − truth).
    blanked = 0
    if blank_keys and level in blank_keys:
        contaminant = blank_keys[level] - t_keys
        blanked = len(f_keys & contaminant)
        f_keys = f_keys - contaminant
    tp, fp, fn = len(f_keys & t_keys), len(f_keys - t_keys), len(t_keys - f_keys)
    return {"experiment": experiment, "sample": sample, "level": level, "mode": "set",
            "n_truth": len(t_keys), "n_found": len(f_keys), "TP": tp, "FP": fp, "FN": fn,
            "blank_subtracted": blanked,
            "fdr_pct": calculate_fdr(f_keys, t_keys), "tpr_pct": calculate_tpr(f_keys, t_keys)}


def _score_protein_level(truth, obs, experiment, sample, q_value_max, exclude, rule: str) -> dict:
    t = _scorable_truth(truth, experiment, sample, exclude)
    o = _scorable_obs(obs, experiment, sample, exclude, q_value_max).dropna(subset=["protein_group"])
    truth_prot = {tok for f in t["protein_names"] for tok in _tokens(f)}

    tp_groups = fp_groups = 0
    found_members: set[str] = set()
    for _, rows in o.groupby("protein_group"):
        members = {tok for f in rows["protein_names"] for tok in _tokens(f)}
        if not members:
            continue
        found_members |= members
        hit = bool(members & truth_prot)
        is_tp = (members <= truth_prot) if rule == "strict" else hit
        tp_groups += int(is_tp)
        fp_groups += int(not is_tp)

    n_groups = tp_groups + fp_groups
    detected = truth_prot & found_members  # recall counts truth proteins, no member-expansion
    # NOTE: protein-row precision (TP/FP/n_found/fdr_pct) is in GROUP units; recall
    # (n_truth_proteins/detected_proteins/tpr_pct) is in PROTEIN units. They are
    # deliberately different bases and do not form one confusion matrix.
    return {"experiment": experiment, "sample": sample, "level": "protein",
            "mode": f"group-credit:{rule}",
            "n_truth": len(truth_prot), "n_found": n_groups,
            "TP": tp_groups, "FP": fp_groups, "FN": len(truth_prot) - len(detected),
            "n_truth_proteins": len(truth_prot), "detected_proteins": len(detected),
            "fdr_pct": round(100.0 * fp_groups / n_groups, 3) if n_groups else 0.0,
            "tpr_pct": round(100.0 * len(detected) / len(truth_prot), 3) if truth_prot else 0.0}


def score_fdr(truth: pd.DataFrame, obs: pd.DataFrame,
              set_levels=("ion", "peptide"), protein_rules=("any", "strict"),
              q_value_max: float = 0.01, exclude_scopes=EXCLUDE_SCOPES,
              blank_keys=None) -> pd.DataFrame:
    """FDR/TPR table per ``(experiment, sample, level, mode)``.

    Scored independently per ``(experiment, sample)`` — truth is never unioned
    across experiments or across A/B. ``blank_keys`` (from :func:`blank_keysets`)
    subtracts IDs found in the blank noise source from the FP side.
    """
    pairs = sorted(set(zip(truth["experiment"], truth["sample"]))
                   & set(zip(obs["experiment"], obs["sample"])))
    rows = []
    for experiment, sample in pairs:
        for level in set_levels:
            rows.append(_score_set_level(truth, obs, experiment, sample, level,
                                         q_value_max, exclude_scopes, blank_keys))
        for rule in protein_rules:
            rows.append(_score_protein_level(truth, obs, experiment, sample,
                                             q_value_max, exclude_scopes, rule))
    return pd.DataFrame(rows)


def plot_fdr_tpr(scores: pd.DataFrame, out: Path,
                 title: str = "PlasmaBENCH P5 — simulated-spike-in FDR / TPR") -> None:
    """Grouped bars: empirical FDR% and recall% per level, one colour per sample."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11})
    df = scores.copy()
    df["label"] = df["level"] + df["mode"].map(
        lambda m: "" if m == "set" else f"\n({m.split(':')[-1]})")
    labels = list(dict.fromkeys(df["label"]))
    samples = sorted(df["sample"].unique())
    palette = {"A": "#4C72B0", "B": "#DD8452"}
    x = np.arange(len(labels))
    w = 0.8 / max(1, len(samples))

    fig, (axf, axt) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    for ax, metric, name in ((axf, "fdr_pct", "empirical FDR %"), (axt, "tpr_pct", "recall (TPR) %")):
        for i, s in enumerate(samples):
            d = df[df["sample"] == s].set_index("label").reindex(labels)
            ax.bar(x + (i - (len(samples) - 1) / 2) * w, d[metric], w,
                   label=f"Sample {s}", color=palette.get(s, None))
        ax.set_xticks(x); ax.set_xticklabels(labels)
        ax.set_ylabel(name); ax.set_title(name)
        ax.grid(axis="y", ls=":", alpha=0.4)
    axf.axhline(1.0, color="red", ls="--", lw=1, label="1% target")
    axf.legend(); axt.legend()
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


__all__ = ["calculate_fdr", "calculate_tpr", "score_fdr", "plot_fdr_tpr", "blank_keysets"]
