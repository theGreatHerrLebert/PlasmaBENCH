"""Translate search-engine peptide notations into TimSim's UniMod sequence.

TimSim's `from_findings` loader feeds each sequence to the Rust PeptideSequence
parser, which accepts only standard residues with inline `[UNIMOD:n]` tags
(plus an optional leading N-term tag) and *aborts the process* on anything
else. So we must (a) translate each engine's notation faithfully and
(b) drop rows that still don't validate, rather than let them reach Rust.

Two engines are supported:

  - **DIA-NN** — `Modified.Sequence` carries inline `(UniMod:n)` tokens
    (N-term mods already lead the sequence). Translation is a notation swap:
    `(UniMod:n)` → `[UNIMOD:n]`.

  - **FragPipe** — the stripped `Peptide` lives in one column and the mods in
    a separate `Assigned Modifications` column (`"1M(15.9949), N-term(42.0106)"`).
    This mirrors the proven parser in
    timsim-bench/timsim_bench/benchmark/helpers.py: mass→UniMod via sagepy's
    `mass_to_mod`, tag emitted *after* the residue, N-term acetyl pulled to the
    front as `[UNIMOD:1]`.
"""
from __future__ import annotations

import re
from typing import List

import numpy as np

# Standard residues + inline [UNIMOD:n], with an optional leading N-term tag.
# Mirrors what the Rust PeptideSequence parser accepts in practice (the
# from_findings path used by MHCbench emits leading N-term tags this way).
_VALID_SEQUENCE_RE = re.compile(
    r"^(\[UNIMOD:\d+\])?([ACDEFGHIKLMNPQRSTVWY](\[UNIMOD:\d+\])?)+$"
)

_DIANN_TOKEN = re.compile(r"\(UniMod:(\d+)\)", re.IGNORECASE)


def diann_modseq_to_unimod(modified_sequence: str) -> str:
    """DIA-NN `Modified.Sequence` → TimSim UniMod sequence.

    Examples:
        >>> diann_modseq_to_unimod("AAC(UniMod:4)DEFK")
        'AAC[UNIMOD:4]DEFK'
        >>> diann_modseq_to_unimod("(UniMod:1)SDKPDM(UniMod:35)AEIEK")
        '[UNIMOD:1]SDKPDM[UNIMOD:35]AEIEK'
    """
    if not isinstance(modified_sequence, str) or not modified_sequence:
        return modified_sequence
    # DIA-NN occasionally wraps sequences in underscores — strip them.
    seq = modified_sequence.strip().strip("_")
    return _DIANN_TOKEN.sub(lambda m: f"[UNIMOD:{m.group(1)}]", seq)


# --- FragPipe: Peptide + Assigned Modifications -----------------------------
# Carried over from timsim-bench/benchmark/helpers.py so the two stay in sync.

def _parse_pattern(mod: str) -> List[str]:
    """Match a FragPipe assigned-mod token like `1M(15.9949)`."""
    match = re.match(r"^(\d+)([A-Za-z])\(([\d.]+)\)$", mod)
    return list(match.groups()) if match else []


def _parse_mods(sequence: str, mods: str):
    """Return (index, replacement) edits for a FragPipe Assigned Modifications cell."""
    r_list = []
    if not isinstance(mods, str) or not mods:
        return r_list
    for mod in mods.split(", "):
        if mod.startswith("N-term(42"):
            # N-terminal acetylation → lead with [UNIMOD:1] on residue 0.
            r_list.append((0, f"[UNIMOD:1]{sequence[0]}"))
        else:
            parts = _parse_pattern(mod)
            if parts:
                index, aa, mass = parts
                tag = mass_to_mod(int(np.round(float(mass))))
                r_list.append((int(index) - 1, aa + tag))
    return r_list


def fragpipe_to_unimod(sequence: str, assigned_modifications: str) -> str:
    """FragPipe stripped `Peptide` + `Assigned Modifications` → UniMod sequence.

    Examples:
        >>> fragpipe_to_unimod("AACDEFK", "3C(57.0215)")
        'AAC[UNIMOD:4]DEFK'
        >>> fragpipe_to_unimod("SDKPDMAEIEK", "6M(15.9949), N-term(42.0106)")
        '[UNIMOD:1]SDKPDM[UNIMOD:35]AEIEK'
    """
    residues = {i: aa for i, aa in enumerate(sequence)}
    for index, replacement in _parse_mods(sequence, assigned_modifications):
        residues[index] = replacement
    return "".join(residues[i] for i in range(len(sequence)))


def is_valid_timsim_sequence(sequence: str) -> bool:
    """True if `sequence` will parse in TimSim's PeptideSequence (won't abort)."""
    return bool(isinstance(sequence, str) and _VALID_SEQUENCE_RE.match(sequence))


# `mass_to_mod` is imported lazily so importing this module for the DIA-NN path
# (or for `is_valid_timsim_sequence`) does not require the heavy sagepy stack.
def mass_to_mod(mass: int) -> str:  # noqa: D103
    from sagepy.core.peptide import mass_to_mod as _m2m
    return _m2m(mass)


__all__ = [
    "diann_modseq_to_unimod",
    "fragpipe_to_unimod",
    "is_valid_timsim_sequence",
]
