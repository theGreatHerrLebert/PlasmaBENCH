# Stage 1 results

Scored metrics land here once `make stage1-eval` is implemented:
- `metrics.json` — species-resolved recall, true FDR, and log2 B/A ratio-recovery
  error per tool.
- per-tool tables under `<tool>/`.

The B/A ratio plot (log2 ratio vs. abundance, coloured by species, with the
expected lines at 0 / +1.585 / −1.0) is the headline figure — the LFQbench-style
view of how tightly each pipeline recovers the known spike-in ratios.
