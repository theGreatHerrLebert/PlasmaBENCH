#!/usr/bin/env bash
# Regenerate every quant-dependent report figure from the current results/ tree.
# Quantity basis = raw Precursor.Quantity (loader default; see
# docs/finding-diann-1.8-vs-2.5-quant.md). FDR figures are q-value based and are
# NOT regenerated here (unaffected by the quantity-column choice).
#
# Usage:  source .venv/bin/activate && bash scripts/make_report_figures.sh
set -euo pipefail
cd "$(dirname "$0")/.."

FIG=results/figures
S1=simulations/stage1/dia
S2=simulations/stage2/dia-fix
R18_S1=results/diann-1.8-AB-v3/report.tsv
R25_S1=results/diann-2.5-AB-v3/report.parquet
R18_S2=results/diann-1.8-s2fix-r080/report.tsv
R25_S2=results/diann-2.5-s2fix-r080/report.parquet
R18_S2_081=results/diann-1.8-s2fix-r081/report.tsv
R25_S2_081=results/diann-2.5-s2fix-r081/report.parquet

echo "### §2  Stage-1 ratio panels"
python scripts/plot_stage1_ratios.py --report "$R25_S1" --sim-dir "$S1" \
  --experiment stage1-diann-2.5 --level protein --protein-quant raw_sum --out "$FIG/stage1_diann25_ratios_protein.png"
python scripts/plot_stage1_ratios.py --report "$R25_S1" --sim-dir "$S1" \
  --experiment stage1-diann-2.5 --level ion --out "$FIG/stage1_diann25_ratios_ion.png"
python scripts/plot_stage1_ratios.py --report "$R18_S1" --sim-dir "$S1" \
  --experiment stage1-diann-1.8 --level protein --protein-quant raw_sum --out "$FIG/stage1_diann18_ratios_protein.png"
python scripts/plot_stage1_ratios.py --report "$R18_S1" --sim-dir "$S1" \
  --experiment stage1-diann-1.8 --level ion --out "$FIG/stage1_diann18_ratios_ion.png"

echo "### §3  Stage-1 quant correlation"
python scripts/plot_stage1_quant.py --report "$R25_S1" --sim-dir "$S1" \
  --experiment stage1-diann-2.5 --out "$FIG/stage1_diann25_quant.png"
python scripts/plot_stage1_quant.py --report "$R18_S1" --sim-dir "$S1" \
  --experiment stage1-diann-1.8 --out "$FIG/stage1_diann18_quant.png"

echo "### §6  Stage-1 true-abundance + ratio-bias (1.8 vs 2.5)"
python scripts/plot_true_abundance.py --report-a "$R18_S1" --label-a "DIA-NN 1.8" \
  --report-b "$R25_S1" --label-b "DIA-NN 2.5" --sim-dir "$S1" --level ion \
  --out "$FIG/stage1_trueabund_ion_18_vs_25.png"
python scripts/plot_ratio_bias.py --report-a "$R18_S1" --label-a "DIA-NN 1.8" \
  --report-b "$R25_S1" --label-b "DIA-NN 2.5" --sim-dir "$S1" --level ion \
  --out "$FIG/stage1_ratiobias_ion_18_vs_25.png"

echo "### §7  Stage-2 (fixed) ratio panels + ratio-bias + cross-stage abundance"
python scripts/plot_stage1_ratios.py --report "$R18_S2" --sim-dir "$S2" \
  --experiment stage2-diann-1.8 --level protein --protein-quant raw_sum --out "$FIG/stage2fix_18_ratios_protein.png"
python scripts/plot_stage1_ratios.py --report "$R25_S2" --sim-dir "$S2" \
  --experiment stage2-diann-2.5 --level protein --protein-quant raw_sum --out "$FIG/stage2fix_25_ratios_protein.png"
# Second background (r081) — no figure embedded, but reproduces the §7b "stable across
# 080/081" claim (both engines, raw_sum protein, anchor ≈ 0 on both). Outputs to /tmp.
echo "  -- r081 (stability check for §7b table) --"
python scripts/plot_stage1_ratios.py --report "$R18_S2_081" --sim-dir "$S2" \
  --experiment stage2-diann-1.8-r081 --level protein --protein-quant raw_sum --out /tmp/s2fix_18_r081.png
python scripts/plot_stage1_ratios.py --report "$R25_S2_081" --sim-dir "$S2" \
  --experiment stage2-diann-2.5-r081 --level protein --protein-quant raw_sum --out /tmp/s2fix_25_r081.png
python scripts/plot_ratio_bias.py --report-a "$R18_S2" --label-a "DIA-NN 1.8" \
  --report-b "$R25_S2" --label-b "DIA-NN 2.5" --sim-dir "$S2" --level ion \
  --out "$FIG/stage2fix_18v25_ratiobias.png"
python scripts/plot_true_abundance.py --report-a "$R18_S1" --label-a "Stage 1 (blank)" \
  --report-b "$R18_S2" --label-b "Stage 2 (real plasma)" --sim-dir "$S1" --sim-dir-b "$S2" \
  --level ion --out "$FIG/stage1_vs_stage2fix_trueabund.png"

echo "### §8  1.8-vs-2.5 quant comparison (raw headline + normalization mechanism)"
python scripts/plot_quant_18v25.py

echo "### §9  real-vs-SIM compression capstone (needs results/diann-1.8-realG + SIM stage reports)"
python scripts/plot_real_vs_sim_compression.py || \
  echo "  (skipped: real-data search results/diann-1.8-realG not present)"
python scripts/plot_lowabund_underquant.py || \
python scripts/plot_real_ratios.py --report results/diann-1.8-realG/report.tsv --label "DIA-NN 1.8" --out "$FIG/real_ratios_diann18.png" || echo "  (skipped: real 1.8 report missing)"
python scripts/plot_real_ratios.py --report results/diann-2.5-realG/report.parquet --label "DIA-NN 2.5" --out "$FIG/real_ratios_diann25.png" || echo "  (skipped: real 2.5 report missing)"
  echo "  (skipped §9b: SIM Stage-2 report/blueprint not present)"
echo "### §9a real-data cross-engine ID overlap (1.8 vs 2.x)"
python scripts/plot_id_overlap.py --report-a results/diann-1.8-realG/report.tsv --label-a "DIA-NN 1.8" \
  --report-b results/diann-2.5-realG/report.parquet --label-b "DIA-NN 2.5" \
  --title "Real PYE1 ID overlap — 1.8 vs 2.5 (q<0.01, 12 runs pooled)" \
  --out "$FIG/real_overlap_18_vs_25.png" || echo "  (skipped: real reports missing)"
python scripts/plot_id_overlap.py --report-a results/diann-1.8-realG/report.tsv --label-a "DIA-NN 1.8" \
  --report-b results/diann-2.6-realG/report.parquet --label-b "DIA-NN 2.6" \
  --title "Real PYE1 ID overlap — 1.8 vs 2.6 (q<0.01, 12 runs pooled)" \
  --out "$FIG/real_overlap_18_vs_26.png" || echo "  (skipped: real reports missing)"

echo "DONE — figures in $FIG"
