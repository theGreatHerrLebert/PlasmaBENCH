# DIA-NN — config notes

_TBD._ Add the `.cfg` / command line used, a `VERSION` file, and any parameters
that could not be set to the unified value.

Seed relevance: `build_seed_from_report.py --engine diann` reads `report.tsv` /
`report.parquet`. Columns consumed: `Protein.Group`/`Protein.Ids`,
`Modified.Sequence` (`(UniMod:n)` notation), `Precursor.Charge`,
`Precursor.Quantity`, `RT` (minutes → converted to seconds), `IM` (1/K0).
Confirm these names against the DIA-NN version Ute used.
