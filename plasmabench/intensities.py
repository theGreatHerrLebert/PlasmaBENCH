"""Synthetic per-peptide baseline intensities for the `imposed` ratio mode.

The `imposed` seed reproduces what TimSim's HYE / proteome_mix pipeline does for
intensity: it does NOT use the search engine's reported (and ratio-compressed)
intensity. Instead each peptide draws a baseline abundance from the Tenzer
rank-abundance curve, and the A/B fold change comes purely from a per-species
dilution factor (see plasmabench/species.py PYE1_COMPOSITION). The result is a
clean, uncompressed ground truth.

Ported verbatim from rustims imspy_simulation/timsim/jobs/simulate_proteins.py
(`get_tenzer_hokey` + `assign_events`) so the two stay in sync:

    events = upscale_factor * tenzer_curve[random_rank]

The one change for seeding is determinism: the random rank is drawn from a
per-peptide RNG keyed on (sequence, charge) via a stable hash, so an independent
Sample-A and Sample-B invocation assign the SAME baseline to a shared peptide —
which is what makes the imposed fold change exact (mirrors HYE's
`is_fold_change_condition`, where condition B reuses condition A's baseline).
"""
from __future__ import annotations

import zlib

import numpy as np

_CURVE_LEN = int(1e4)


def get_tenzer_curve() -> np.ndarray:
    """The Tenzer rank-abundance curve (length 1e4), highest abundance at rank 1.

    Verbatim from rustims `get_tenzer_hokey()`.
    """
    delta_exp = 0.005
    delta = 0.0011
    exp_red = 600
    target_start = 1e4

    rank = np.linspace(1, _CURVE_LEN, _CURVE_LEN)
    exp_values = 1 / np.exp(rank / exp_red)
    cumulative_exp_sum = np.cumsum(exp_values)
    return target_start * (2 ** (-delta * rank)) * (2 ** (-delta_exp * cumulative_exp_sum))


def _stable_seed(sequence: str, charge) -> int:
    """Deterministic 32-bit seed for a peptide ion, stable across processes.

    Python's built-in hash() is salted per-process, so it would give Sample A
    and Sample B different baselines — use crc32 instead.
    """
    return zlib.crc32(f"{sequence}|{charge}".encode("utf-8")) & 0xFFFFFFFF


def baseline_intensity(sequence: str, charge, *, upscale_factor: float,
                       curve: np.ndarray | None = None) -> float:
    """Per-peptide baseline abundance = upscale_factor * tenzer[random_rank].

    The rank is drawn deterministically from (sequence, charge) so the same ion
    gets the same baseline in every sample — the basis of an exact fold change.
    """
    if curve is None:
        curve = get_tenzer_curve()
    rng = np.random.RandomState(_stable_seed(sequence, charge))
    # rank drawn from [1, _CURVE_LEN) to match rustims `assign_events`
    # (np.random.uniform(1, 1e4)) exactly — index 0 (the single highest-abundance
    # point) is intentionally not sampled, same as upstream. Does not affect A/B
    # ratios (baseline is shared between samples).
    rank = int(rng.uniform(1, _CURVE_LEN))
    rank = min(rank, _CURVE_LEN - 1)
    return float(upscale_factor) * float(curve[rank])


__all__ = ["get_tenzer_curve", "baseline_intensity"]
