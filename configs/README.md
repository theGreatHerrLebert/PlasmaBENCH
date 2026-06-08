# Search-engine configs

One subdirectory per tool. Each must contain:
- the search config in the tool's native format
- a `VERSION` file with the exact tool version used
- a short `NOTES.md` calling out parameters that could not be set to the
  unified value

Tools (start with the two already in Ute's report set):
- `dia-nn/`
- `fragpipe/`
- (later) `spectronaut/`, `maxquant/`

Settings must be unified across tools where the knob exists (precursor/fragment
mass tolerance, peptide length range, variable + fixed modifications, FDR
target, decoy export). Plasma samples are alkylated, so carbamidomethyl-C is a
fixed modification (unlike the MHCbench immunopeptidome configs).
