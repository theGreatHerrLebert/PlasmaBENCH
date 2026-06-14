"""Mechanism panel: the SIM low-abundance over-separation is asymmetric under-quantification
of the LOWER-intensity member of each A/B pair (report §9b).

Decomposes the over-separation into per-member quant recovery log2(observed/truth), binned by
TRUE abundance (blueprint), per species. The higher member is recovered accurately at all
abundances; the lower member is under-quantified at low abundance — and the higher-minus-lower
gap equals the over-separation. SIM Stage 2 (DIA-NN 1.8). The global truth->observed scale
cancels in the higher-vs-lower difference, so the curves are directly interpretable.
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
from plasmabench.plots.loader import build_manifest, load_truth, load_observations

FIG = Path("results/figures/sim_lowabund_underquant.png")
REP = "results/diann-1.8-s2fix-r080/report.tsv"
SIM = "simulations/stage2/dia-fix"
NBINS = 8
SPECIES = ["ECOLI", "YEAST"]


def _wide(df, col):
    df = df.copy(); df["feat"] = df.sequence_modified + df.charge.astype(str)
    df = df[df[col] > 0]
    w = df.groupby(["feat", "sample"])[col].mean().unstack()
    sp = df.groupby("feat")["species"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None)
    w["species"] = w.index.map(sp)
    return w


def main():
    FIG.parent.mkdir(parents=True, exist_ok=True)
    m = build_manifest(REP, SIM, "s2", "background_unknown")
    truth = pd.concat([load_truth(r, m) for r in m.runs], ignore_index=True)
    obs = load_observations(REP, m)
    # q<0.01 filter on observations before joining (score only confident IDs).
    for qc in ("q_value", "lib_q_value", "pg_q_value", "lib_pg_q_value"):
        if qc in obs.columns:
            obs = obs[obs[qc].isna() | (obs[qc] <= 0.01)]
    j = _wide(truth, "truth_intensity").join(_wide(obs, "observed_intensity"),
                                             lsuffix="_t", rsuffix="_o", how="inner")
    j = j.dropna(subset=["A_t", "B_t", "A_o", "B_o"])
    j = j[j.species_t.isin(SPECIES)].copy()
    for s in ("A", "B"):
        j[f"rec_{s}"] = np.log2(j[f"{s}_o"]) - np.log2(j[f"{s}_t"])
    ctr = np.median(np.r_[j.rec_A, j.rec_B]); j.rec_A -= ctr; j.rec_B -= ctr
    hi_is_B = j.B_t >= j.A_t
    j["rec_hi"] = np.where(hi_is_B, j.rec_B, j.rec_A)
    j["rec_lo"] = np.where(hi_is_B, j.rec_A, j.rec_B)
    j["ab"] = 0.5 * np.log2(j.A_t * j.B_t)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharey=True)
    for ax, sp in zip(axes, SPECIES):
        d = j[j.species_t == sp].copy(); d["oct"] = pd.qcut(d.ab, NBINS, labels=False, duplicates="drop")
        x = np.arange(1, NBINS + 1)
        hi = d.groupby("oct")["rec_hi"].median(); lo = d.groupby("oct")["rec_lo"].median()
        ax.plot(x, hi, "-o", color="#4C72B0", ms=5, lw=1.8, label="higher member (big peak)")
        ax.plot(x, lo, "-s", color="#C44E52", ms=5, lw=1.8, label="lower member (small peak)")
        ax.fill_between(x, lo, hi, color="0.8", alpha=0.5, label="= over-separation")
        ax.axhline(0, color="0.5", lw=1.0, ls=(0, (4, 3)))
        ax.set_title(sp.title(), fontsize=11)
        ax.set_xlabel("TRUE-abundance octile  (1 = lowest → 8 = highest)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("quant recovery  log2(observed / truth)")
    axes[0].legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle("Cause of the SIM over-separation: the LOW-intensity peak is under-quantified at "
                 "low abundance; the big peak is accurate (SIM Stage 2, DIA-NN 1.8)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"wrote {FIG}")


if __name__ == "__main__":
    main()
