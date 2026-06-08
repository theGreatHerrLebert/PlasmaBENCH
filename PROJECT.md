# PlasmaBENCH — kanban

Mirrors the GitHub Project board so the plan is reviewable in the repo.

## Backlog
- Collect GitHub usernames (distleu, tenzer); add as repo collaborators.
- Decide repo visibility (start **private**) once the team agrees.
- Decide raw-data hosting (fileshare link vs. DVC/S3 mirror).
- CI: lint configs, validate the EVIDENT manifest.

## Stage 1 — fully simulated PYE1 (in flight)
- [ ] Receive Sample A / Sample B DIA-NN **and** FragPipe reports from Ute.
- [ ] Receive ≥1 TimsTOF DIA-PASEF **blank** `.d` (same method as the study); ideally a plasma-only `.d` too.
- [ ] Confirm PYE1 ratios (A = 90/2/8, B = 90/6/4 — plasma/yeast/E. coli).
- [ ] Validate `build_seed_from_report.py` against Ute's actual report columns/engine versions.
- [ ] Build `seed_sampleA.csv` / `seed_sampleB.csv`; sanity-check species split and intensity ranges.
- [ ] Render configs and run TimSim DIA for both samples.
- [ ] Document the exact `timsim` invocation + rustims commit per run.

## Stage 2 — SIM spiked onto real plasma (planned)
- [ ] Acquire / locate a real plasma-only `.d` matching the study method.
- [ ] Wire `superimpose_on_reference` against it; compare to the experimental TimsTOF runs.

## Software setup (cross-cutting)
- [ ] Unify search parameters across DIA-NN / FragPipe (ppm, length, mods, FDR); flag untunable knobs per tool.
- [ ] Confirm decoy export per tool (entrapment-style validation).
- [ ] Capture license/version per tool in `configs/<tool>/VERSION`.

## Decisions logged
- 2026-06-08 — Repo scaffolded from the MHCbench pattern. TimsTOF only (no Orbitrap arm). Seed route = `from_findings` from real reports; ratios baked in, no `proteome_mix`.
