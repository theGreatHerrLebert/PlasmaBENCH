"""PlasmaBENCH analysis/plotting layer.

LFQbench-style panels (ported from the PYE multicenter paper's R scripts) plus
ground-truth-only panels (FDR, quant correlation, controlled-jitter precision)
that the paper structurally cannot produce — TimSim emits a per-simulation
blueprint of exact truth.

Everything hangs on two normalized tables built by `loader.py`:

- ``truth``        — from each run's ``synthetic_data.db`` blueprint
- ``observations`` — from the engine (DIA-NN/FragPipe) report

resolved against a shared run :class:`~plasmabench.plots.loader.RunManifest`.
See ``plan.md`` for the full design and the decisions behind it.
"""
