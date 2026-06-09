#!/usr/bin/env bash
# Build the PlasmaBENCH search FASTA: human (plasma background) + S. cerevisiae
# EC1118 + E. coli K12 — the three proteomes behind Ute's PYE1 reports, picked to
# match the organism tags the seeds carry (_HUMAN / _YEAS8 / _ECOLI).
#
# Verified from the report accessions (scripts identify them via UniProt):
#   human   taxId 9606    reviewed Swiss-Prot   (P02768=ALBU_HUMAN, ...)
#   ecoli   taxId 83333   reviewed Swiss-Prot   K12, proteome UP000000625
#   yeast   taxId 643680  TrEMBL                 EC1118, proteome UP000000286 (_YEAS8)
#
# Output (gitignored): data/raw/fasta/plasma_yeast_ecoli.fasta
#
# NOTE: this is a best-effort reconstruction of Ute's search database. Prefer her
# EXACT FASTA if she can share it; run scripts/qc_fasta_coverage.py afterwards to
# measure how many seed peptides this database actually covers.
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR="$REPO_ROOT/data/raw/fasta"
mkdir -p "$OUT_DIR"
STREAM="https://rest.uniprot.org/uniprotkb/stream"

PROV="$OUT_DIR/FASTA_PROVENANCE.txt"
: > "$PROV"

fetch() {  # <label> <query> <outfile> <min_expected>
  local label="$1" query="$2" out="$3" min="$4"
  echo "==> $label"
  local hdr; hdr=$(mktemp)
  curl -sS --fail --retry 3 -G "$STREAM" -D "$hdr" \
    --data-urlencode "query=$query" \
    --data-urlencode "format=fasta" \
    --data-urlencode "includeIsoform=false" \
    -o "$out"
  local n; n=$(grep -c '^>' "$out")
  # Fail loudly on an HTTP-200-but-empty/short response — otherwise a UniProt
  # outage or query change would silently produce a degraded benchmark DB.
  if (( n < min )); then
    echo "ERROR: '$label' returned $n sequences (< expected $min). Query: $query" >&2
    rm -f "$hdr"; exit 3
  fi
  local rel; rel=$(grep -i '^x-uniprot-release' "$hdr" | tr -d '\r' || true)
  echo "    $n sequences -> $out   [$rel]"
  printf '%s\t%s seqs\t%s\tquery=%s\n' "$label" "$n" "${rel:-release?}" "$query" >> "$PROV"
  rm -f "$hdr"
}

fetch "human (reviewed Swiss-Prot, 9606)" \
      "(organism_id:9606) AND (reviewed:true)" "$OUT_DIR/_human.fasta" 18000
fetch "E. coli K12 (reviewed, UP000000625)" \
      "(proteome:UP000000625) AND (reviewed:true)" "$OUT_DIR/_ecoli.fasta" 4000
# Use the FULL set of EC1118 (_YEAS8) entries by organism, NOT the redundancy-
# reduced reference proteome UP000000286 (4514): the full set is ~5984 and the
# extra ~1470 entries are needed to cover the seed yeast peptides (the reference
# proteome alone only covers ~77%). See scripts/qc_fasta_coverage.py.
fetch "S. cerevisiae EC1118 (organism 643680, all _YEAS8)" \
      "(organism_id:643680)" "$OUT_DIR/_yeast.fasta" 5000

OUT="$OUT_DIR/plasma_yeast_ecoli.fasta"
cat "$OUT_DIR/_human.fasta" "$OUT_DIR/_yeast.fasta" "$OUT_DIR/_ecoli.fasta" > "$OUT"
rm -f "$OUT_DIR/_human.fasta" "$OUT_DIR/_ecoli.fasta" "$OUT_DIR/_yeast.fasta"

echo
echo "Combined: $OUT"
echo "  total sequences: $(grep -c '^>' "$OUT")"
echo "  by tag: HUMAN=$(grep -c '_HUMAN' "$OUT")  YEAS8=$(grep -c '_YEAS8' "$OUT")  ECOLI=$(grep -c '_ECOLI' "$OUT")"
