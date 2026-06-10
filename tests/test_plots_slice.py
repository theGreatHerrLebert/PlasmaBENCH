"""Integration tests for the plasmabench.plots loader + P6 quant panel.

Runs against the REAL Stage-1 fixtures (the Sample-A/B blueprints and the
DIA-NN 2.5 report) — no pytest required:

    python tests/test_plots_slice.py

Skips (exit 0) if the local fixtures are absent, so it is safe in CI without data.
"""
from __future__ import annotations

import sqlite3
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plasmabench.plots.loader import (  # noqa: E402
    RunSpec, RunManifest, load_truth, load_observations, TRUTH_COLUMNS, OBS_COLUMNS,
)
from plasmabench.plots.quant import build_quant_table, correlations, JOIN_KEYS  # noqa: E402
from plasmabench.plots import ratios as R  # noqa: E402
from plasmabench.plots import fdr as F  # noqa: E402

REPORT = ROOT / "results/diann-2.5-AB-v3/report.parquet"
SIMDIR = ROOT / "simulations/stage1/dia"
RUN_A = "PLB-S1-DIA-sampleA-imposed"
RUN_B = "PLB-S1-DIA-sampleB-imposed"
DB_A = SIMDIR / f"sampleA_imposed/{RUN_A}/synthetic_data.db"
DB_B = SIMDIR / f"sampleB_imposed/{RUN_B}/synthetic_data.db"

PASSED = 0


def check(cond: bool, msg: str) -> None:
    global PASSED
    if not cond:
        raise AssertionError(msg)
    PASSED += 1
    print(f"  ok: {msg}")


def manifest() -> RunManifest:
    return RunManifest(runs=[
        RunSpec("stage1", "A", RUN_A, DB_A),
        RunSpec("stage1", "B", RUN_B, DB_B),
    ])


def test_manifest_rejects_duplicates() -> None:
    print("test_manifest_rejects_duplicates")
    try:
        RunManifest(runs=[RunSpec("e", "A", RUN_A, DB_A), RunSpec("e", "B", RUN_A, DB_B)])
        raise SystemExit("FAIL: duplicate run_id was accepted")
    except ValueError as e:
        check("duplicate run_id" in str(e), "duplicate run_id raises ValueError")
    m = manifest()
    try:
        m.resolve("not-a-run")
        raise SystemExit("FAIL: unknown run resolved")
    except KeyError:
        check(True, "unknown run_id raises KeyError")
    check(m.scope_for("HUMAN") == "simulated", "human scope=simulated by default")
    check(m.scope_for("YEAST") == "simulated", "yeast scope=simulated")
    check(m.scope_for(None) == "unscorable", "ambiguous species → unscorable")
    # YE-only condition: human becomes background_unknown
    m2 = RunManifest(runs=[RunSpec("e", "A", RUN_A, DB_A)], human_scope="background_unknown")
    check(m2.scope_for("HUMAN") == "background_unknown", "YE-only: human=background_unknown")
    check(m2.scope_for("ECOLI") == "simulated", "YE-only: ecoli still simulated")


def test_truth_contract() -> None:
    print("test_truth_contract")
    t = load_truth(RunSpec("stage1", "A", RUN_A, DB_A), manifest())
    check(list(t.columns) == TRUTH_COLUMNS, "truth has exactly the contract columns")
    check(t.duplicated(JOIN_KEYS).sum() == 0, "truth precursor key is unique")
    # transmitted is a strict subset (untransmitted precursors exist).
    n_tx = int(t["transmitted"].sum())
    check(0 < n_tx < len(t), f"transmitted is a strict subset ({n_tx}/{len(t)})")
    check((t["truth_intensity"] > 0).all(), "all truth_intensity > 0")
    check(t["species"].isin(["HUMAN", "YEAST", "ECOLI"]).any(), "species assigned")
    check(t["sequence_modified"].str.contains(r"\[UNIMOD:", regex=True).any(),
          "blueprint sequences carry UniMod tags")
    return t


def test_observations_contract() -> None:
    print("test_observations_contract")
    o = load_observations(REPORT, manifest())
    check(set(o.columns) <= set(OBS_COLUMNS), "observations columns ⊆ contract")
    check(o.duplicated(JOIN_KEYS).sum() == 0, "observation precursor key is unique")
    check(set(o["run_id"].unique()) == {RUN_A, RUN_B}, "both runs resolved via manifest")
    check(o["species"].isin(["HUMAN", "YEAST", "ECOLI"]).any(), "species assigned from Protein.Names")
    # unknown run must fail loudly
    bad = RunManifest(runs=[RunSpec("stage1", "A", RUN_A, DB_A)])  # B missing
    try:
        load_observations(REPORT, bad)
        raise SystemExit("FAIL: report with unmapped run loaded")
    except KeyError:
        check(True, "report run absent from manifest raises KeyError")
    return o


def test_quant_join_and_corr(t: pd.DataFrame, o: pd.DataFrame) -> None:
    print("test_quant_join_and_corr")
    tB = load_truth(RunSpec("stage1", "B", RUN_B, DB_B), manifest())
    truth = pd.concat([t, tB], ignore_index=True)
    merged = build_quant_table(truth, o, q_value_max=0.01)
    check(len(merged) > 40_000, f"join matched many precursors ({len(merged):,})")
    check(merged.duplicated(JOIN_KEYS).sum() == 0, "join stayed one-to-one")
    check(merged["truth_intensity"].gt(0).all() and merged["observed_intensity"].gt(0).all(),
          "matched points are strictly positive")
    c = correlations(merged)
    check((c["spearman"] > 0.7).all(), "every (sample,species) Spearman > 0.7 (imposed seeds)")
    check(set(c["species"]) == {"HUMAN", "YEAST", "ECOLI"}, "all three species scored")


def test_ratio_convention_single_api() -> None:
    print("test_ratio_convention_single_api")
    # the one display transform, exercised on data and on expected lines together
    import numpy as np
    check(R.A_OVER_B.to_display(0.5) == -0.5, "A/B flips sign of canonical log2(B/A)")
    check(R.B_OVER_A.to_display(0.5) == 0.5, "B/A is identity on canonical")
    # expected lines derive from the SAME transform → paper axis values
    check(abs(R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA["YEAST"]) - (-np.log2(3))) < 1e-9,
          "A/B yeast expected == -log2(3)")
    check(abs(R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA["ECOLI"]) - 1.0) < 1e-9,
          "A/B ecoli expected == +1")
    check(R.A_OVER_B.to_display(R.EXPECTED_LOG2_BA["HUMAN"]) == 0.0, "A/B human expected == 0")
    # B/A is the mirror
    check(abs(R.B_OVER_A.to_display(R.EXPECTED_LOG2_BA["YEAST"]) - np.log2(3)) < 1e-9,
          "B/A yeast expected == +log2(3)")


def test_ratio_table_recovers_oracle(o: pd.DataFrame) -> None:
    print("test_ratio_table_recovers_oracle")
    wide, anchor = R.build_ratio_table(o, level="protein", aggregation="linear_mean",
                                       human_anchor=True)
    check(len(wide) > 3000, f"protein ratio table built ({len(wide):,} proteins)")
    check(set(wide["species"]) == {"HUMAN", "YEAST", "ECOLI"}, "all three species present")
    summ = R.ratio_summary(wide, R.A_OVER_B).set_index("species")
    # observed median within 0.4 log2 of the oracle (imposed seeds, slight compression)
    for sp in ("YEAST", "ECOLI"):
        d = abs(summ.loc[sp, "observed_log2"] - summ.loc[sp, "expected_log2"])
        check(d < 0.4, f"{sp} observed within 0.4 log2 of expected ({d:.2f})")
    check(abs(summ.loc["HUMAN", "observed_log2"]) < 1e-6, "human anchored to 0")

    # all four levels build, with feature counts nesting protein < peptide <= ion
    n = {}
    for lvl in ("protein", "peptide", "modified_peptide", "ion"):
        w, _ = R.build_ratio_table(o, level=lvl)
        n[lvl] = len(w)
        s = R.ratio_summary(w, R.A_OVER_B).set_index("species")
        check(abs(s.loc["YEAST", "observed_log2"] - s.loc["YEAST", "expected_log2"]) < 0.4,
              f"{lvl}: yeast within 0.4 log2 of expected")
    check(n["protein"] < n["peptide"] <= n["ion"],
          f"feature counts nest: protein {n['protein']} < peptide {n['peptide']} <= ion {n['ion']}")
    check(n["modified_peptide"] <= n["ion"], "modified_peptide collapses charges vs ion")
    # 'precursor' is an alias for 'ion'
    check(len(R.build_ratio_table(o, level="precursor")[0]) == n["ion"],
          "precursor == ion alias")

    # protein_quant="raw_sum" (raw Σ-precursor) also recovers the oracle on the fixture
    wsum, _ = R.build_ratio_table(o, level="protein", protein_quant="raw_sum")
    ssum = R.ratio_summary(wsum, R.A_OVER_B).set_index("species")
    for sp in ("YEAST", "ECOLI"):
        d = abs(ssum.loc[sp, "observed_log2"] - ssum.loc[sp, "expected_log2"])
        check(d < 0.4, f"raw_sum {sp} observed within 0.4 log2 of expected ({d:.2f})")


def test_fdr_math() -> None:
    print("test_fdr_math")
    found, truth = {1, 2, 3}, {2, 3, 4}
    check(abs(F.calculate_fdr(found, truth) - 100 / 3) < 1e-2, "FDR = FP/(FP+TP)")
    check(abs(F.calculate_tpr(found, truth) - 200 / 3) < 1e-2, "TPR = TP/(TP+FN)")
    check(F.calculate_fdr(set(), {1}) == 0.0, "empty found → 0% FDR")


def test_fdr_scoring(o: pd.DataFrame) -> None:
    print("test_fdr_scoring")
    truth = pd.concat([
        load_truth(RunSpec("stage1", "A", RUN_A, DB_A), manifest()),
        load_truth(RunSpec("stage1", "B", RUN_B, DB_B), manifest()),
    ], ignore_index=True)
    s = F.score_fdr(truth, o, q_value_max=0.01)
    check(set(s["sample"]) == {"A", "B"}, "both samples scored, A/B not unioned")
    ion = s[(s["level"] == "ion") & (s["sample"] == "A")].iloc[0]
    # transmitted-only denominator: ion truth count == transmitted truth in sample A
    n_tx_A = int(truth[(truth["sample"] == "A") & truth["transmitted"]
                       & truth["truth_scope"].eq("simulated") & truth["species"].notna()].shape[0])
    check(ion["n_truth"] == n_tx_A, f"ion truth set = transmitted simulated truth ({n_tx_A})")
    check((s["fdr_pct"] < 5).all(), "empirical FDR < 5% everywhere (q<0.01)")
    check((s["tpr_pct"] > 50).all(), "recall > 50% everywhere")
    # strict group-credit is never more lenient than any-member
    for smp in ("A", "B"):
        anyf = s[(s["level"] == "protein") & (s["mode"] == "group-credit:any") & (s["sample"] == smp)]["fdr_pct"].iloc[0]
        strf = s[(s["level"] == "protein") & (s["mode"] == "group-credit:strict") & (s["sample"] == smp)]["fdr_pct"].iloc[0]
        check(strf >= anyf, f"sample {smp}: strict FDR ≥ any-member FDR")
    # YE-only condition: human dropped from scoring (background_unknown)
    ye_manifest = RunManifest(runs=[RunSpec("stage1", "A", RUN_A, DB_A),
                                    RunSpec("stage1", "B", RUN_B, DB_B)],
                              human_scope="background_unknown")
    tru_ye = pd.concat([load_truth(r, ye_manifest) for r in ye_manifest.runs], ignore_index=True)
    obs_ye = load_observations(REPORT, ye_manifest)
    s_ye = F.score_fdr(tru_ye, obs_ye, q_value_max=0.01)
    ion_ye = s_ye[(s_ye["level"] == "ion") & (s_ye["sample"] == "A")].iloc[0]
    check(ion_ye["n_truth"] < ion["n_truth"], "YE-only: human excluded shrinks truth set")


# ---------------------------------------------------------------------------
# Synthetic defect-coverage tests (no fixtures needed — always run).
# Each targets a specific bug found in the Codex code review.
# ---------------------------------------------------------------------------
from plasmabench.plots.loader import _collapse_observation_duplicates, OBS_COLUMNS  # noqa: E402


def _make_blueprint(path: str) -> None:
    """Tiny synthetic_data.db exercising loader truth-collapse edge cases."""
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE peptides (peptide_id INT, protein_id INT, sequence TEXT, "
                "protein TEXT, decoy INT, events REAL)")
    con.execute("CREATE TABLE ions (ion_id INT, peptide_id INT, charge INT, relative_abundance REAL)")
    con.execute("CREATE TABLE fragment_ions (ion_id INT)")
    # PEPTIDEK/2 reachable via two same-species proteins; one ion transmitted (100),
    # one NOT (200). SAMPLERK/3 untransmitted-only. DECOYK is a decoy (excluded).
    con.executemany("INSERT INTO peptides VALUES (?,?,?,?,?,?)", [
        (1, 10, "PEPTIDEK", "AAA_YEAS8", 0, 100.0),
        (2, 11, "PEPTIDEK", "BBB_YEAS8", 0, 200.0),
        (3, 12, "SAMPLERK", "CCC_ECOLI", 0, 50.0),
        (4, 13, "DECOYK",   "DDD_ECOLI", 1, 999.0),
    ])
    con.executemany("INSERT INTO ions VALUES (?,?,?,?)", [
        (1, 1, 2, 1.0), (2, 2, 2, 1.0), (3, 3, 3, 1.0), (4, 4, 2, 1.0)])
    con.executemany("INSERT INTO fragment_ions VALUES (?)", [(1,)])  # only ion 1 transmitted
    con.commit(); con.close()


def test_loader_truth_collapse_synthetic() -> None:
    print("test_loader_truth_collapse_synthetic")
    with tempfile.TemporaryDirectory() as d:
        db = Path(d) / "synthetic_data.db"
        _make_blueprint(str(db))
        m = RunManifest(runs=[RunSpec("e", "A", "r", db)])
        t = load_truth(RunSpec("e", "A", "r", db), m).set_index("sequence_modified")
        check("DECOYK" not in t.index, "decoy excluded from truth")
        pep = t.loc["PEPTIDEK"]
        check(bool(pep["transmitted"]), "PEPTIDEK transmitted (one ion in fragment_ions)")
        # transmitted-aware: only the transmitted ion's 100 counts, NOT 100+200
        check(abs(pep["truth_intensity"] - 100.0) < 1e-9,
              f"truth_intensity = transmitted-only (100, not 300); got {pep['truth_intensity']}")
        check("AAA_YEAS8" in pep["protein_names"] and "BBB_YEAS8" in pep["protein_names"],
              "protein_names UNIONs both same-species proteins (P5 recall)")
        check(pep["species"] == "YEAST", "union resolves to YEAST")
        check(pd.notna(pep["sim_protein_id"]), "sim_protein_id populated (not NaN)")
        check(not bool(t.loc["SAMPLERK"]["transmitted"]), "untransmitted-only precursor flagged")


def test_obs_collapse_nullkey_and_dedup() -> None:
    print("test_obs_collapse_nullkey_and_dedup")
    base = dict(experiment="e", sample="A", run_id="r", protein_names="X_HUMAN",
                species="HUMAN", truth_scope="simulated", protein_group="P", software="DIA-NN")
    # duplicate (seq,charge) → median of observed_intensity
    df = pd.DataFrame([
        {**base, "sequence_modified": "AAA", "charge": 2, "observed_intensity": 10.0, "q_value": 0.001},
        {**base, "sequence_modified": "AAA", "charge": 2, "observed_intensity": 30.0, "q_value": 0.002},
        {**base, "sequence_modified": "BBB", "charge": 2, "observed_intensity": 5.0, "q_value": 0.001},
    ])
    collapsed = _collapse_observation_duplicates(df)
    check(len(collapsed) == 2, "duplicate (seq,charge) collapsed to one row")
    aaa = collapsed[collapsed.sequence_modified == "AAA"].iloc[0]
    check(abs(aaa["observed_intensity"] - 20.0) < 1e-9, "duplicate intensity → median (20)")
    # a null join key must raise, not silently vanish
    bad = pd.concat([df, pd.DataFrame([{**base, "sequence_modified": "AAA", "charge": 2,
                                        "observed_intensity": 1.0, "q_value": None,
                                        "run_id": None}])], ignore_index=True)
    try:
        _collapse_observation_duplicates(bad)
        raise SystemExit("FAIL: null join key silently accepted")
    except ValueError:
        check(True, "null join key raises instead of dropping rows")


def _fdr_frame(rows, truth=True):
    cols_t = ["experiment", "sample", "run_id", "sequence_modified", "charge", "sequence",
              "species", "truth_scope", "transmitted", "protein_names"]
    cols_o = ["experiment", "sample", "run_id", "sequence_modified", "charge", "sequence",
              "species", "truth_scope", "protein_names", "protein_group", "q_value"]
    return pd.DataFrame(rows, columns=cols_t if truth else cols_o)


def test_fdr_experiment_scope_and_ambiguous_fp() -> None:
    print("test_fdr_experiment_scope_and_ambiguous_fp")
    # Two experiments share sample 'A'. K2 is truth only in E2; E1 finding K2 must be FP.
    truth = _fdr_frame([
        ["E1", "A", "r1", "K1", 2, "K1", "YEAST", "simulated", True, "Y_YEAS8"],
        ["E2", "A", "r2", "K2", 2, "K2", "YEAST", "simulated", True, "Y_YEAS8"],
    ], truth=True)
    obs = _fdr_frame([
        ["E1", "A", "r1", "K1", 2, "K1", "YEAST", "simulated", "Y_YEAS8", "Y", 0.001],   # TP
        ["E1", "A", "r1", "K2", 2, "K2", "YEAST", "simulated", "Y_YEAS8", "Y", 0.001],   # FP (truth only in E2)
        ["E1", "A", "r1", "K3", 2, "K3", None,    "unscorable", "Z;W",     "Z", 0.001],  # ambiguous FP, must count
        ["E1", "A", "r1", "KH", 2, "KH", "HUMAN", "background_unknown", "H_HUMAN", "H", 0.001],  # excluded
    ], truth=False)
    s = F.score_fdr(truth, obs, set_levels=("ion",), protein_rules=())
    e1 = s[(s["experiment"] == "E1") & (s["level"] == "ion")].iloc[0]
    check(e1["TP"] == 1, "K1 is the only TP in E1 (no cross-experiment credit)")
    check(e1["FP"] == 2, "K2 (other-experiment truth) AND K3 (ambiguous) both count as FP")
    check(e1["n_found"] == 3, "background_unknown (KH) excluded from found, others kept")


def test_render_key_collision_and_proteome_guard() -> None:
    print("test_render_key_collision_and_proteome_guard")
    from scripts.render_configs import render
    base = ('reference_path  = "/blank.d"\n'
            'reference_in_memory  = false\n'
            'reference_noise_intensity_max = 150000\n'
            'add_real_data_noise = true\n'
            '# superimpose_on_reference = false\n'
            'findings_path = "x"\nsave_path = "y"\nexperiment_name = "z"\n'
            'proteome_mix = true\n')
    out = render(base, sample="A", seed_csv="S", save_path="SP", experiment_name="EN",
                 reference_path="/plasma.d", superimpose=True)
    check('reference_path  = "/plasma.d"' in out, "reference_path repointed")
    check("reference_in_memory  = false" in out, "reference_in_memory NOT clobbered by prefix match")
    check("reference_noise_intensity_max = 150000" in out, "reference_noise_* NOT clobbered")
    check(out.count('reference_path  =') == 1, "exactly one reference_path line (no duplicate)")
    check("add_real_data_noise           = false" in out, "superimpose flips add_real_data_noise")
    check("superimpose_on_reference      = true" in out, "superimpose uncommented & set true")
    check("proteome_mix = false" in out and "proteome_mix = true" not in out,
          "superimpose ENFORCES proteome_mix=false")


def test_fdr_blank_subtraction() -> None:
    print("test_fdr_blank_subtraction")
    truth = _fdr_frame([
        ["E1", "A", "r1", "K1", 2, "K1", "YEAST", "simulated", True, "Y_YEAS8"],
    ], truth=True)
    obs = _fdr_frame([
        ["E1", "A", "r1", "K1", 2, "K1", "YEAST", "simulated", "Y_YEAS8", "Y", 0.001],  # TP
        ["E1", "A", "r1", "KB", 2, "KB", "YEAST", "simulated", "Y_YEAS8", "Y", 0.001],  # FP, but in blank
        ["E1", "A", "r1", "KF", 2, "KF", "YEAST", "simulated", "Y_YEAS8", "Y", 0.001],  # FP, genuine
    ], truth=False)
    no_blank = F.score_fdr(truth, obs, set_levels=("ion",), protein_rules=())
    fp0 = no_blank[no_blank["level"] == "ion"].iloc[0]["FP"]
    check(fp0 == 2, "without blank: KB and KF both FP")
    # KB is found in the blank → should drop off the FP side; KF stays
    blank = {"ion": {("KB", 2)}}
    with_blank = F.score_fdr(truth, obs, set_levels=("ion",), protein_rules=(), blank_keys=blank)
    row = with_blank[with_blank["level"] == "ion"].iloc[0]
    check(row["FP"] == 1, "blank-explained FP (KB) subtracted; genuine FP (KF) remains")
    check(row["TP"] == 1 and row["blank_subtracted"] == 1, "TP unchanged, 1 blank ID subtracted")


def _abund_frames(experiments=("E",)):
    """Matched truth/obs ion-level frames; ECOLI 'EEEEEEK' is transmitted only in A."""
    trows, orows = [], []
    ec = ["AAAAAAK", "CCCCCCK", "DDDDDDK", "EEEEEEK"]
    for exp in experiments:
        for i, seq in enumerate(ec):
            for samp, inten in (("A", 1000 * (i + 1)), ("B", 500 * (i + 1))):
                tx = not (seq == "EEEEEEK" and samp == "B")  # E4 transmitted only in A
                trows.append(dict(experiment=exp, sample=samp, sequence_modified=seq, charge=2,
                                  sequence=seq, species="ECOLI", truth_scope="simulated",
                                  transmitted=tx, truth_intensity=float(inten)))
                orows.append(dict(experiment=exp, sample=samp, run_id=f"{exp}{samp}",
                                  sequence_modified=seq, charge=2,
                                  sequence=seq, species="ECOLI", truth_scope="simulated",
                                  observed_intensity=float(inten), q_value=0.001,
                                  protein_names="X_ECOLI", protein_group="X"))
        for seq in ["HHHHHHK", "IIIIIIK"]:  # human, for the anchor
            for samp in ("A", "B"):
                trows.append(dict(experiment=exp, sample=samp, sequence_modified=seq, charge=2,
                                  sequence=seq, species="HUMAN", truth_scope="simulated",
                                  transmitted=True, truth_intensity=2000.0))
                orows.append(dict(experiment=exp, sample=samp, run_id=f"{exp}{samp}",
                                  sequence_modified=seq, charge=2,
                                  sequence=seq, species="HUMAN", truth_scope="simulated",
                                  observed_intensity=2000.0, q_value=0.001,
                                  protein_names="H_HUMAN", protein_group="H"))
    return pd.DataFrame(trows), pd.DataFrame(orows)


def test_abundance_denominator_and_guard() -> None:
    print("test_abundance_denominator_and_guard")
    from plasmabench.plots.abundance import build_true_abundance
    truth, obs = _abund_frames()
    t = build_true_abundance(truth, obs, level="ion", n_bins=2)
    ec = t[t["species"] == "ECOLI"]
    check(int(ec["n_truth"].sum()) == 4,
          "A-only-transmitted feature stays in the denominator (transmitted-in-≥1, not both)")
    check((ec["sensitivity"] == 1.0).all(), "all eligible+detected ecoli → sensitivity 1.0")
    # empty-but-schema'd: a species with no rows must not break downstream
    from plasmabench.plots.abundance import ABUND_COLUMNS
    check(list(t.columns) == ABUND_COLUMNS, "table has the declared columns even with 1 species")
    # multi-experiment must be rejected (no cross-experiment feature credit)
    t2, o2 = _abund_frames(("E1", "E2"))
    try:
        build_true_abundance(t2, o2, level="ion", n_bins=2)
        raise SystemExit("FAIL: multi-experiment input not rejected")
    except ValueError:
        check(True, "multi-experiment input raises")


def test_tokens_no_pipe_split() -> None:
    print("test_tokens_no_pipe_split")
    check(F._tokens("sp|P12345|ALBU_HUMAN") == ["sp|P12345|ALBU_HUMAN"],
          "_tokens does NOT shatter a sp|acc|name header on '|'")
    check(F._tokens("A_YEAS8;B_YEAS8") == ["A_YEAS8", "B_YEAS8"], "_tokens splits on ';'")


def test_protein_quant_raw_sum() -> None:
    print("test_protein_quant_raw_sum")
    # One yeast protein group, two precursors, A→B is 3× by raw intensity. PG.MaxLFQ is
    # set INCONSISTENT with the raw sum (equal A/B) so the two modes give different answers:
    # raw_sum must read the summed precursor intensities, maxlfq must read PG.MaxLFQ.
    def rows(sample, inten):
        return [dict(experiment="e", sample=sample, run_id=f"r{sample}",
                     protein_names="ALPHA_YEAS8", protein_group="G1", pg_maxlfq=500.0,
                     species="YEAST", sequence_modified=sm, charge=2,
                     observed_intensity=inten, q_value=0.001)
                for sm in ("PEPA", "PEPB")]
    obs = pd.DataFrame(rows("A", 100.0) + rows("B", 300.0))
    w_raw, _ = R.build_ratio_table(obs, level="protein", protein_quant="raw_sum",
                                   human_anchor=False)
    w_mlfq, _ = R.build_ratio_table(obs, level="protein", protein_quant="maxlfq",
                                    human_anchor=False)
    # raw_sum: log2(sum_B / sum_A) = log2(600/200) = +log2 3
    check(abs(float(w_raw["log2_ba"].iloc[0]) - np.log2(3)) < 1e-9,
          f"raw_sum sums precursors → log2(B/A)=+log2 3 ({float(w_raw['log2_ba'].iloc[0]):+.3f})")
    # maxlfq: equal PG.MaxLFQ A/B → 0  (proves the modes read different columns)
    check(abs(float(w_mlfq["log2_ba"].iloc[0])) < 1e-9,
          "maxlfq reads PG.MaxLFQ (equal A/B → 0), not the raw sum")
    try:
        R.build_ratio_table(obs, level="protein", protein_quant="bogus")
        check(False, "invalid protein_quant should raise ValueError")
    except ValueError:
        check(True, "invalid protein_quant raises ValueError")


def run_synthetic() -> None:
    test_loader_truth_collapse_synthetic()
    test_obs_collapse_nullkey_and_dedup()
    test_fdr_experiment_scope_and_ambiguous_fp()
    test_fdr_blank_subtraction()
    test_abundance_denominator_and_guard()
    test_render_key_collision_and_proteome_guard()
    test_tokens_no_pipe_split()
    test_protein_quant_raw_sum()


def test_abundance_curves(o: pd.DataFrame) -> None:
    print("test_abundance_curves")
    from plasmabench.plots.abundance import build_true_abundance
    from plasmabench.plots import ratios as RR
    truth = pd.concat([load_truth(RunSpec("stage1", "A", RUN_A, DB_A), manifest()),
                       load_truth(RunSpec("stage1", "B", RUN_B, DB_B), manifest())],
                      ignore_index=True)
    t = build_true_abundance(truth, o, level="ion")
    check(set(t["species"]) == {"HUMAN", "YEAST", "ECOLI"}, "true-abundance table has all species")
    check(t["sensitivity"].between(0, 1).all(), "sensitivity in [0,1]")
    # sensitivity rises with true abundance (LOD shape): high decile >> low decile for ecoli
    ec = t[t["species"] == "ECOLI"].sort_values("x")
    check(ec.iloc[-1]["sensitivity"] > ec.iloc[0]["sensitivity"] + 0.5,
          "ecoli sensitivity increases strongly with true abundance (LOD curve)")
    # bias near expected at the top decile
    top = ec.iloc[-1]["median_ratio"]
    check(abs(top - RR.EXPECTED_LOG2_BA["ECOLI"] * RR.A_OVER_B.sign) < 0.3,
          f"ecoli high-abundance bias near expected ({top:+.2f})")
    # engine-axis bias helper also builds
    wide, _ = RR.build_ratio_table(o, level="ion")
    b = RR.bias_by_abundance(wide, RR.A_OVER_B)
    check(set(b["species"]) == {"HUMAN", "YEAST", "ECOLI"} and (b["n"] > 0).all(),
          "bias_by_abundance builds per-species deciles")


def main() -> int:
    run_synthetic()  # fixture-free, always run
    for f in (REPORT, DB_A, DB_B):
        if not f.exists():
            print(f"SKIP: fixture-dependent tests (fixture missing: {f})")
            print(f"\nSYNTHETIC PASSED ({PASSED} checks)")
            return 0
    test_manifest_rejects_duplicates()
    t = test_truth_contract()
    o = test_observations_contract()
    test_quant_join_and_corr(t, o)
    test_ratio_convention_single_api()
    test_ratio_table_recovers_oracle(o)
    test_abundance_curves(o)
    test_fdr_math()
    test_fdr_scoring(o)
    print(f"\nALL PASSED ({PASSED} checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
