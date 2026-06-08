"""Generic table / sequence utilities.

Cherry-picked from
  /scratch/timsim-demo/SUBMISSION/timsim-bench/timsim_bench/benchmark/utility.py
keeping only helpers that are not phospho- or plot-specific, so PlasmaBENCH's
seed/analysis scripts can reuse them without depending on timsim-bench.
"""
from __future__ import annotations

import re
from typing import Iterable

import numpy as np


def extract_protein_id(protein_string: str) -> str | None:
    """Pull the UniProt accession out of a `sp|P12345|FOO_HUMAN` style string."""
    match = re.search(r"\|([^|]+)\|", protein_string)
    return match.group(1) if match else None


def replace_I_with_L(sequence: str) -> str:
    """Replace I with L outside `[UNIMOD:...]` brackets (I/L-equivalence)."""
    def replace_match(match: re.Match[str]) -> str:
        if match.group(1):
            return match.group(1)
        return match.group(2).replace("I", "L")

    modified = re.sub(r"(\[.*?\])|([^\[\]]+)", replace_match, sequence)
    return modified.replace("UNLMOD", "UNIMOD")


def calculate_fdr(found: Iterable, target: Iterable) -> float:
    """Empirical FDR % between a found set and a known-truth set."""
    found_set = found if isinstance(found, set) else set(found)
    target_set = target if isinstance(target, set) else set(target)
    fp = len(found_set - target_set)
    tp = len(found_set & target_set)
    total = fp + tp
    return 0.0 if total == 0 else float(np.round((fp / total) * 100, 2))


__all__ = ["extract_protein_id", "replace_I_with_L", "calculate_fdr"]
