"""Organism assignment and the PYE1 expected-ratio table.

The three-species design lives or dies on a clean organism map: every
quantified peptide must be attributable to human (plasma), yeast, or
E. coli, and the few that aren't (shared/ambiguous, contaminants) must be
excluded from the ratio oracle rather than silently miscounted.

Assignment is by UniProt entry-name suffix (`_HUMAN`, `_YEAS8`, `_ECOLI`),
which is what DIA-NN and FragPipe carry through in their `Protein.Names`
column when the search FASTA is the standard SwissProt download per
organism. If Ute's FASTA uses different organism tags, extend `SPECIES_TAGS`.

NOTE: tags live in `Protein.Names` (e.g. `ALBU_HUMAN`), NOT `Protein.Group`
(bare accessions like `P02768`). Feed the names column here.

The yeast proteome Ute used is the *strain-specific* S. cerevisiae download,
whose UniProt OS code is `YEAS8` (strain ATCC 204508 / S288c) — so the tag is
`_YEAS8`, not the canonical `_YEAST`. Both are mapped to YEAST so a mixed
FASTA still resolves. Contaminant entries (`_CONTA`, the cRAP set) match no
tag and therefore return None — correctly excluded from the ratio oracle.
"""
from __future__ import annotations

import re

# UniProt OS tags → canonical species label. E. coli SwissProt entries are
# `*_ECOLI`; human `*_HUMAN`; yeast appears as `*_YEAS8` (strain-specific) or
# `*_YEAST` (canonical) — both map to YEAST.
SPECIES_TAGS: dict[str, str] = {
    "_HUMAN": "HUMAN",
    "_YEAS8": "YEAST",
    "_YEAST": "YEAST",
    "_ECOLI": "ECOLI",
}

# PYE1 nominal composition (fraction of total protein mass per sample).
PYE1_COMPOSITION: dict[str, dict[str, float]] = {
    "A": {"HUMAN": 0.90, "YEAST": 0.02, "ECOLI": 0.08},
    "B": {"HUMAN": 0.90, "YEAST": 0.06, "ECOLI": 0.04},
}

# Expected Sample-B / Sample-A intensity ratio per species — the quant oracle.
EXPECTED_RATIO_B_OVER_A: dict[str, float] = {
    sp: PYE1_COMPOSITION["B"][sp] / PYE1_COMPOSITION["A"][sp]
    for sp in PYE1_COMPOSITION["A"]
}  # HUMAN 1.0, YEAST 3.0, ECOLI 0.5


def assign_species(protein_field: str) -> str | None:
    """Return 'HUMAN' | 'YEAST' | 'ECOLI', or None if absent/ambiguous.

    A peptide mapping to more than one organism (protein groups that span
    species) is ambiguous and returns None — such peptides must not enter
    the ratio oracle.

    Examples:
        >>> assign_species("sp|P02768|ALBU_HUMAN")
        'HUMAN'
        >>> assign_species("P00330_YEAST;P00331_YEAST")
        'YEAST'
        >>> assign_species("sp|P02768|ALBU_HUMAN;sp|P0A6F5|CH60_ECOLI") is None
        True
    """
    if not isinstance(protein_field, str) or not protein_field:
        return None
    found: set[str] = set()
    for tag, label in SPECIES_TAGS.items():
        if tag in protein_field:
            found.add(label)
    if len(found) == 1:
        return next(iter(found))
    return None  # zero matches or cross-species → ambiguous


def expected_log2_ratio(species: str) -> float | None:
    """log2(Sample B / Sample A) expected for a species, or None if unknown."""
    import math
    r = EXPECTED_RATIO_B_OVER_A.get(species)
    return None if r is None else math.log2(r)


__all__ = [
    "SPECIES_TAGS",
    "PYE1_COMPOSITION",
    "EXPECTED_RATIO_B_OVER_A",
    "assign_species",
    "expected_log2_ratio",
]
