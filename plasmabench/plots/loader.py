"""The two-table loader: ``truth`` (blueprint) and ``observations`` (engine report).

Every panel in :mod:`plasmabench.plots` is a pure function of these two tables,
joined on ``(experiment, sample, run_id, sequence_modified, charge)``. The design
and the review that shaped it live in ``plan.md``; the load-bearing invariants:

- **A run manifest is the join authority.** The blueprint (`synthetic_data.db`) is
  one DB per run with no run/sample/experiment column, and DIA-NN ``Run`` stems do
  not match ``.d``/dir names. So identity comes from an explicit manifest, which
  fails on any missing or duplicate mapping — never from path inference.
- **`transmitted` truth is the primary identification universe.** Not every
  simulated ion reaches ``fragment_ions``; untransmitted precursors are unfindable
  and must not count as false negatives. ``transmitted = ion ∈ fragment_ions``.
- **`truth_intensity = events × relative_abundance`** — ``events`` carries the
  baked-in A/B dilution; ``relative_abundance`` splits a peptide across charges.
- **`peptides.protein_id` is an internal integer** → kept as ``sim_protein_id``.
  The joinable ``protein_id`` is the normalized protein token (entry name, e.g.
  ``CFAB_HUMAN``, which is what our `from_findings` blueprint stores and what
  DIA-NN ``Protein.Names`` carries); a full ``sp|acc|name`` header collapses to its
  accession. (Protein-level joins are P5; the precursor panels do not need it.)
- **Species** via :func:`plasmabench.species.assign_species` on the protein-names
  field (handles ``_YEAS8`` and multi-species→None).
"""
from __future__ import annotations

import re
import sqlite3
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plasmabench.species import assign_species  # noqa: E402
from plasmabench.peptides import diann_modseq_to_unimod  # noqa: E402

# Canonical column contracts (see plan.md).
TRUTH_COLUMNS = [
    "experiment", "sample", "run_id", "sequence", "sequence_modified", "charge",
    "precursor_id", "sim_protein_id", "protein_id", "protein_names", "species",
    "truth_intensity", "transmitted", "truth_scope",
]
OBS_COLUMNS = [
    "experiment", "sample", "run_id", "software", "sequence", "sequence_modified",
    "charge", "precursor_id", "protein_id", "protein_names", "protein_group",
    "species", "truth_scope", "observed_intensity", "pg_maxlfq",
    "q_value", "lib_q_value", "pg_q_value", "lib_pg_q_value",
]

_SPECIES_SCOPE = {"YEAST": "simulated", "ECOLI": "simulated"}  # human filled per-manifest


@dataclass(frozen=True)
class RunSpec:
    """One physical run: a `.d` + its `synthetic_data.db`, with canonical identity."""
    experiment: str
    sample: str          # 'A' | 'B'
    run_id: str          # MUST equal the engine report's `Run` value
    blueprint_db: Path


@dataclass
class RunManifest:
    """The single source of truth mapping ``run_id`` → canonical identity + blueprint.

    Fails loudly on duplicate ``run_id`` (the ambiguity Codex flagged). ``human_scope``
    is ``"simulated"`` for full-SIM runs and ``"background_unknown"`` for the
    YE-only-SIM-over-real-plasma condition (human has no blueprint truth there).
    """
    runs: list[RunSpec]
    human_scope: str = "simulated"
    _by_run: dict[str, RunSpec] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.human_scope not in ("simulated", "background_unknown"):
            raise ValueError(f"human_scope must be simulated|background_unknown, got {self.human_scope!r}")
        by_run: dict[str, RunSpec] = {}
        for r in self.runs:
            if r.run_id in by_run:
                raise ValueError(f"duplicate run_id in manifest: {r.run_id!r}")
            by_run[r.run_id] = r
        self._by_run = by_run

    def resolve(self, run_id: str) -> RunSpec:
        try:
            return self._by_run[run_id]
        except KeyError:
            raise KeyError(
                f"run_id {run_id!r} not in manifest (known: {sorted(self._by_run)})"
            ) from None

    def scope_for(self, species: str | None) -> str:
        if species == "HUMAN":
            return self.human_scope
        return _SPECIES_SCOPE.get(species, "unscorable")


def _strip_unimod(seq: str) -> str:
    """`[UNIMOD:4]` tags removed → bare residue string."""
    return re.sub(r"\[UNIMOD:\d+\]", "", seq)


def _union(series) -> str:
    """`;`-joined unique protein tokens, order-preserving (matches build_seed)."""
    return ";".join(dict.fromkeys(str(x) for x in series if pd.notna(x) and str(x)))


def _canonical_protein(protein_field: str) -> str:
    """Normalize a protein token to its joinable id.

    `from_findings` blueprints and DIA-NN `Protein.Names` store entry names
    (`CFAB_HUMAN`); a full `sp|P02768|ALBU_HUMAN` header collapses to its accession.
    Multi-protein groups keep the first token (species ambiguity is handled
    separately via `assign_species` on the full field).
    """
    if not isinstance(protein_field, str) or not protein_field:
        return ""
    first = re.split(r"[;,]", protein_field.strip())[0].strip()
    parts = first.split("|")
    return parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else first


def load_truth(run: RunSpec, manifest: RunManifest, exclude_decoys: bool = True) -> pd.DataFrame:
    """Build the ``truth`` table for one run from its blueprint ``synthetic_data.db``.

    Precursor key is ``(sequence_modified, charge)``; duplicates (a modified
    sequence reachable via >1 internal peptide_id) are collapsed by **summing
    truth_intensity** and OR-ing ``transmitted``, with species kept only if
    unambiguous across the collapsed rows.
    """
    db = Path(run.blueprint_db)
    if not db.exists():
        raise FileNotFoundError(f"blueprint not found for run {run.run_id!r}: {db}")
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        peptides = pd.read_sql_query(
            "SELECT peptide_id, protein_id, sequence, protein, decoy, events FROM peptides", con
        )
        ions = pd.read_sql_query(
            "SELECT ion_id, peptide_id, charge, relative_abundance FROM ions", con
        )
        transmitted_ids = pd.read_sql_query(
            "SELECT DISTINCT ion_id FROM fragment_ions", con
        )["ion_id"]
    finally:
        con.close()

    df = ions.merge(peptides, on="peptide_id", how="left", validate="many_to_one")
    df["transmitted"] = df["ion_id"].isin(set(transmitted_ids))
    if exclude_decoys:
        df = df[df["decoy"] == 0]

    df["sequence_modified"] = df["sequence"]
    df["truth_intensity"] = df["events"].astype(float) * df["relative_abundance"].astype(float)
    df["protein_names"] = df["protein"]
    df["sim_protein_id"] = df["protein_id"]            # internal integer, kept as-is
    # transmitted-only intensity, so untransmitted duplicates can't leak into an
    # observable precursor's truth_intensity (see collapse below).
    df["_ti_tx"] = df["truth_intensity"].where(df["transmitted"], 0.0)

    # Collapse to one row per (sequence_modified, charge): UNION the protein names
    # (a peptide may map to >1 same-species protein — keep them all for P5 recall),
    # intensity = transmitted ions' sum if any transmitted else the full sum.
    out = df.groupby(["sequence_modified", "charge"], sort=False).agg(
        sim_protein_id=("sim_protein_id", "first"),
        protein_names=("protein_names", _union),
        transmitted=("transmitted", "any"),
        _ti_all=("truth_intensity", "sum"),
        _ti_tx=("_ti_tx", "sum"),
    ).reset_index()
    out["truth_intensity"] = np.where(out["transmitted"], out["_ti_tx"], out["_ti_all"])
    # species / canonical id from the UNION (assign_species resolves cross-species → None)
    out["species"] = out["protein_names"].map(assign_species)
    out["protein_id"] = out["protein_names"].map(_canonical_protein)
    out["sequence"] = out["sequence_modified"].map(_strip_unimod)
    out = out.drop(columns=["_ti_all", "_ti_tx"])
    out["experiment"] = run.experiment
    out["sample"] = run.sample
    out["run_id"] = run.run_id
    out["precursor_id"] = out["sequence_modified"] + out["charge"].astype(str)
    out["truth_scope"] = out["species"].map(manifest.scope_for)

    out = out[TRUTH_COLUMNS]
    _assert_unique_precursor(out, "truth")
    return out


def load_observations(report_path: str | Path, manifest: RunManifest,
                      software: str = "DIA-NN", exclude_decoys: bool = True,
                      quant_col: str | None = None) -> pd.DataFrame:
    """Build the ``observations`` table from a DIA-NN report (parquet or tsv).

    Every ``Run`` is resolved through the manifest (fails on an unknown run).
    Precursor key ``(sequence_modified, charge)`` per run is asserted unique.

    ``quant_col`` overrides which DIA-NN column becomes ``observed_intensity``.
    Default is RAW ``Precursor.Quantity`` (the benchmark standard — see the qcol
    comment below and docs/finding-diann-1.8-vs-2.5-quant.md); pass
    ``quant_col="Precursor.Normalised"`` for the engine's normalized quantity, or
    ``"Ms1.Area"`` etc. Used by the quant-method diagnostics to compare bases.
    """
    report_path = Path(report_path)
    if report_path.suffix.lower() in (".parquet", ".pq"):
        df = pd.read_parquet(report_path)
    else:
        df = pd.read_csv(report_path, sep="\t", low_memory=False)

    if exclude_decoys and "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]

    # Resolve identity through the manifest (raises on unknown run_id).
    runs = {rid: manifest.resolve(rid) for rid in df["Run"].unique()}
    df["experiment"] = df["Run"].map(lambda r: runs[r].experiment)
    df["sample"] = df["Run"].map(lambda r: runs[r].sample)
    df["run_id"] = df["Run"]
    df["software"] = software

    df["sequence_modified"] = df["Modified.Sequence"].map(diann_modseq_to_unimod)
    df["sequence"] = df.get("Stripped.Sequence", df["sequence_modified"].map(_strip_unimod))
    df["charge"] = df["Precursor.Charge"].astype(int)
    df["protein_names"] = df.get("Protein.Names")
    df["protein_group"] = df.get("Protein.Group")
    df["protein_id"] = df.get("Protein.Ids")
    df["species"] = df["protein_names"].map(assign_species)
    df["truth_scope"] = df["species"].map(manifest.scope_for)

    if quant_col is not None:
        if quant_col not in df.columns:
            raise ValueError(f"quant_col {quant_col!r} not in report columns")
        qcol = quant_col
    else:
        # Benchmark standard = RAW Precursor.Quantity, not Precursor.Normalised.
        # DIA-NN's cross-run normalization assumes most signal is unchanged between
        # runs; on a deliberate spike-in with controlled (simulated) loading that
        # assumption is false, so the .Normalised column injects large per-precursor
        # A/B ratio variance + a compression bias — catastrophically so for 2.5
        # (~7× wider ratio IQR than 1.8). Raw quantity isolates peak extraction from
        # each engine's normalization choice and is the fair cross-engine basis.
        # See docs/finding-diann-1.8-vs-2.5-quant.md. Override via quant_col= for the
        # default-normalized view a user gets out-of-the-box.
        if "Precursor.Quantity" in df.columns:
            qcol = "Precursor.Quantity"
        else:
            qcol = "Precursor.Normalised"
            warnings.warn(
                "Precursor.Quantity absent — falling back to Precursor.Normalised; "
                "ratios will carry the engine's cross-run normalization, NOT the raw "
                "benchmark standard. Results are not directly comparable across engines.",
                RuntimeWarning, stacklevel=2)
    df["observed_intensity"] = df[qcol].astype(float)
    df["pg_maxlfq"] = df.get("PG.MaxLFQ")
    df["q_value"] = df.get("Q.Value")
    df["lib_q_value"] = df.get("Lib.Q.Value")
    df["pg_q_value"] = df.get("PG.Q.Value")
    df["lib_pg_q_value"] = df.get("Lib.PG.Q.Value")
    df["precursor_id"] = df.get("Precursor.Id", df["sequence_modified"] + df["charge"].astype(str))

    out = df[[c for c in OBS_COLUMNS if c in df.columns]].copy()
    out = _collapse_observation_duplicates(out)
    _assert_unique_precursor(out, "observations")
    return out


def read_diann_minimal(path: str | Path, q_value_max: float | None = None) -> pd.DataFrame:
    """Manifest-free DIA-NN report → just the identity + level columns.

    For runs with no blueprint (e.g. the blank noise source): yields
    ``sequence_modified, sequence, charge, protein_group, protein_names, q_value``
    with the UniMod notation bridge applied, decoys removed.
    """
    path = Path(path)
    df = pd.read_parquet(path) if path.suffix.lower() in (".parquet", ".pq") \
        else pd.read_csv(path, sep="\t", low_memory=False)
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    out = pd.DataFrame({
        "sequence_modified": df["Modified.Sequence"].map(diann_modseq_to_unimod),
        "charge": df["Precursor.Charge"].astype(int),
    })
    out["sequence"] = df["Stripped.Sequence"] if "Stripped.Sequence" in df.columns \
        else out["sequence_modified"].map(_strip_unimod)
    out["protein_group"] = df.get("Protein.Group")
    out["protein_names"] = df.get("Protein.Names")
    out["q_value"] = df.get("Q.Value")
    if q_value_max is not None and "q_value" in out.columns:
        out = out[out["q_value"].fillna(1.0) <= q_value_max]
    return out


def _collapse_observation_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Documented duplicate policy: a (run, precursor) seen >1× → median quantity.

    DIA-NN main reports are already one row per (Run, Precursor.Id); this is a
    safety net (e.g. channel rows) rather than the common path.
    """
    keys = ["experiment", "sample", "run_id", "sequence_modified", "charge"]
    null_key = df[keys].isna().any(axis=1)
    if null_key.any():
        raise ValueError(f"observations have {int(null_key.sum())} rows with a null join "
                         f"key {keys} — cannot resolve precursor identity")
    dup = df.duplicated(keys, keep=False)
    if not dup.any():
        return df
    num = {"observed_intensity": "median", "pg_maxlfq": "median",
           "q_value": "min", "lib_q_value": "min", "pg_q_value": "min", "lib_pg_q_value": "min"}
    agg = {c: (num[c] if c in num else "first") for c in df.columns if c not in keys}
    # dropna=False: a null key would have raised above, but never silently drop rows.
    return df.groupby(keys, sort=False, as_index=False, dropna=False).agg(agg)


def _assert_unique_precursor(df: pd.DataFrame, name: str) -> None:
    keys = ["experiment", "sample", "run_id", "sequence_modified", "charge"]
    n_dup = int(df.duplicated(keys).sum())
    if n_dup:
        raise AssertionError(f"{name}: {n_dup} duplicate rows on {keys}")


def _sample_of(run_id: str) -> str:
    """Infer 'A'/'B' from a run id (e.g. ...sampleA...)."""
    m = re.search(r"sample\s*([AB])", run_id, re.IGNORECASE)
    if not m:
        raise ValueError(f"cannot infer sample (A/B) from run_id {run_id!r}")
    return m.group(1).upper()


def _find_blueprint(sim_dir: Path, run_id: str) -> Path:
    hits = list(Path(sim_dir).glob(f"**/{run_id}/synthetic_data.db"))
    if not hits:
        raise FileNotFoundError(f"no synthetic_data.db for run {run_id!r} under {sim_dir}")
    if len(hits) > 1:
        raise ValueError(f"ambiguous blueprint for run {run_id!r}: {hits}")
    return hits[0]


def build_manifest(report_path: str | Path, sim_dir: str | Path, experiment: str,
                   human_scope: str = "simulated") -> RunManifest:
    """Construct a :class:`RunManifest` from an engine report + the SIM output tree.

    Each ``Run`` in the report is mapped to its blueprint at
    ``<sim-dir>/**/<run_id>/synthetic_data.db`` and a sample inferred from the run id.
    """
    report_path = Path(report_path)
    if report_path.suffix.lower() in (".parquet", ".pq"):
        runs_col = pd.read_parquet(report_path, columns=["Run"])["Run"]
    else:
        runs_col = pd.read_csv(report_path, sep="\t", usecols=["Run"])["Run"]
    specs = [
        RunSpec(experiment=experiment, sample=_sample_of(rid), run_id=rid,
                blueprint_db=_find_blueprint(sim_dir, rid))
        for rid in sorted(runs_col.unique())
    ]
    return RunManifest(runs=specs, human_scope=human_scope)


__all__ = [
    "RunSpec", "RunManifest", "load_truth", "load_observations", "build_manifest",
    "TRUTH_COLUMNS", "OBS_COLUMNS",
]
