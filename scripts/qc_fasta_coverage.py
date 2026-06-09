"""QC: does the search FASTA contain the seed peptides?

The seeds are the peptides TimSim will put into the SIM .d. If the FASTA
(searched by DIA-NN/FragPipe) doesn't contain a seed peptide's protein, that
peptide can never be identified — an ARTIFICIAL false negative that has nothing
to do with the software under test. This is the main risk of reconstructing
Ute's database from UniProt rather than using her exact FASTA (strain mismatch,
isoforms, contaminants). Run this to quantify the miss rate per species.

Method: in-silico Trypsin/P digest of the FASTA (missed cleavages 0..N, length
7..30), per species; then check what fraction of each seed's (mod-stripped)
sequences land in the matching species' peptide set.

Usage:
    python scripts/qc_fasta_coverage.py \
        --fasta data/raw/fasta/plasma_yeast_ecoli.fasta \
        --seeds simulations/stage1/seeds/seed_sampleA_imposed.csv \
                simulations/stage1/seeds/seed_sampleB_imposed.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.species import assign_species  # noqa: E402

_UNIMOD = re.compile(r"\[UNIMOD:\d+\]", re.IGNORECASE)


def strip_mods(seq: str) -> str:
    return _UNIMOD.sub("", str(seq))


def digest(seq: str, missed: int, min_len: int, max_len: int):
    """Trypsin/P: cut after every K and R. Yield peptides with 0..missed MC."""
    # cleavage sites = index after each K/R
    cuts = [0] + [i + 1 for i, a in enumerate(seq) if a in "KR"]
    if cuts[-1] != len(seq):
        cuts.append(len(seq))
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(cuts))):
            pep = seq[cuts[i]:cuts[j]]
            if min_len <= len(pep) <= max_len:
                yield pep


def build_peptide_sets(fasta: Path, missed: int, min_len: int, max_len: int):
    sets: dict[str, set] = {"HUMAN": set(), "YEAST": set(), "ECOLI": set()}
    header, sp, buf = None, None, []
    n_prot = 0

    def flush():
        nonlocal n_prot
        if sp in sets and buf:
            n_prot_local = "".join(buf)
            for pep in digest(n_prot_local, missed, min_len, max_len):
                sets[sp].add(pep)

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                sp = assign_species(header)
                buf = []
                n_prot += 1
            else:
                buf.append(line.strip())
    flush()
    return sets, n_prot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--seeds", required=True, nargs="+", type=Path)
    ap.add_argument("--missed-cleavages", type=int, default=1,
                    help="must match the search config (configs/dia-nn/diann-pye1.cfg "
                         "uses 1); a higher value over-counts coverage by treating "
                         "peptides DIA-NN won't search as findable.")
    ap.add_argument("--min-len", type=int, default=7)
    ap.add_argument("--max-len", type=int, default=30)
    args = ap.parse_args()

    print(f"digesting {args.fasta} (Trypsin/P, MC<= {args.missed_cleavages}, "
          f"len {args.min_len}-{args.max_len}) ...")
    sets, n_prot = build_peptide_sets(args.fasta, args.missed_cleavages,
                                      args.min_len, args.max_len)
    print(f"  {n_prot:,} proteins -> peptides per species: "
          + ", ".join(f"{k}={len(v):,}" for k, v in sets.items()))

    seed = pd.concat([pd.read_csv(s) for s in args.seeds], ignore_index=True)
    seed = seed.drop_duplicates(subset=["sequence"])
    seed["bare"] = seed["sequence"].map(strip_mods)
    seed["sp"] = seed["protein"].map(assign_species)

    print(f"\nseed peptides (unique sequences): {len(seed):,}")
    print(f"{'species':7} {'seed pep':>9} {'covered':>9} {'coverage':>9} {'missed':>8}")
    total_cov = total = 0
    for sp in ["HUMAN", "YEAST", "ECOLI"]:
        s = seed[seed["sp"] == sp]
        if not len(s):
            continue
        covered = s["bare"].isin(sets[sp]).sum()
        total_cov += covered
        total += len(s)
        miss = len(s) - covered
        print(f"{sp:7} {len(s):>9,} {covered:>9,} {100*covered/len(s):>8.1f}% {miss:>8,}")
    none = seed["sp"].isna().sum()
    print(f"\noverall coverage: {total_cov:,}/{total:,} ({100*total_cov/total:.1f}%)"
          f"   [{none:,} seed rows had no species tag, ignored]")
    print("\nMisses = seed peptides absent from the reconstructed DB -> would be "
          "artificial false negatives. High yeast miss rate => get Ute's exact "
          "EC1118 FASTA (strain variants).")


if __name__ == "__main__":
    main()
