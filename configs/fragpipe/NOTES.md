# FragPipe — config notes

_TBD._ Add the `.workflow` file, a `VERSION` file, and any parameters that could
not be set to the unified value.

Seed relevance: `build_seed_from_report.py --engine fragpipe` reads the stripped
`Peptide`/`Peptide Sequence` plus a separate `Assigned Modifications` column
(`"6M(15.9949), N-term(42.0106)"`), translated via
`plasmabench/peptides.fragpipe_to_unimod` (mass→UniMod through sagepy's
`mass_to_mod`). Also consumes `Charge`, `Intensity`, `Retention` (seconds),
`Ion Mobility`. Confirm against the FragPipe version Ute used.
