"""Render per-sample TimSim configs from a base template.

The base config uses container-shaped paths:
    /work = repo root
    /data = ${PLB_DATA:-<repo>/data}        (raw blanks + plasma + fastas)
    /sim  = ${PLB_OUT:-<repo>/simulations}  (TimSim outputs)
That works inside the rustims Docker image (volumes mounted by Makefile /
docker-compose). From a host venv (`make venv-source`) those paths don't
exist, so `--path-mode host` rewrites them to host equivalents.

Substitutes SAMPLE_PLACEHOLDER and repoints findings_path / save_path /
experiment_name for each sample.

Usage:
    python scripts/render_configs.py \
        --base simulations/stage1/configs/base-dia.toml \
        --seeds-dir simulations/stage1/seeds \
        --out-dir simulations/stage1/dia \
        --samples A B \
        --path-mode host
"""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


def host_substitutions(repo_root: Path) -> dict[str, str]:
    repo = repo_root.resolve()
    data = os.environ.get("PLB_DATA") or str(repo / "data")
    sim = os.environ.get("PLB_OUT") or str(repo / "simulations")
    return {"/work": str(repo), "/data": data, "/sim": sim}


def rewrite_paths(text: str, subs: dict[str, str]) -> str:
    out = text
    for src, dst in subs.items():
        pattern = re.compile(rf'"({re.escape(src)})((?:/[^"]*)?)"')
        out = pattern.sub(lambda m: f'"{dst}{m.group(2)}"', out)
    return out


def rewrite_value(value: str, subs: dict[str, str]) -> str:
    """Apply container→host prefix subs to a single path value (host mode)."""
    for src, dst in subs.items():
        if value == src or value.startswith(src + "/"):
            return dst + value[len(src):]
    return value


def _toml_key(line: str) -> str | None:
    """The TOML key of a ``key = value`` line (ignoring a leading ``#``), else None.

    Exact-key matching avoids prefix collisions — e.g. ``reference_path`` must not
    also match ``reference_in_memory`` / ``reference_noise_intensity_max``.
    """
    s = line.strip()
    if s.startswith("#"):
        s = s[1:].strip()
    if "=" not in s:
        return None
    return s.split("=", 1)[0].strip()


def render(base_toml: str, *, sample: str, seed_csv: str, save_path: str,
           experiment_name: str, reference_path: str | None = None,
           superimpose: bool = False, gradient_length: str | None = None,
           findings_reference_median: float | None = None,
           intensity_multiplier: float | None = None) -> str:
    """Fill the per-sample fields; optionally switch to the Stage-2 superimpose mode.

    Stage-2 (overlay simulated spike-in onto a real plasma .d) flips three keys:
    ``add_real_data_noise=false``, ``superimpose_on_reference=true`` (the two are
    mutually exclusive), and repoints ``reference_path`` at the plasma run. The
    simulated signal then sits on the *real* human background — so human carries no
    blueprint truth (score with ``--human-scope background_unknown``).

    ``findings_reference_median`` (+ a single ``intensity_multiplier`` shared across
    samples) is the native, preferred way to preserve A/B ratios under from_findings:
    both samples divide events by the SAME reference median instead of their own, so
    median(A)/median(B) no longer distorts the cross-sample ratio (rustims PR #407).
    This replaces the interim per-sample ``--intensity-multiplier`` work-around. Inject
    both as top-level keys right after ``findings_path`` (config.X reads them directly).
    """
    out = base_toml.replace("SAMPLE_PLACEHOLDER", sample)
    new_lines = []
    for line in out.splitlines():
        key = _toml_key(line)
        if key == "findings_path":
            new_lines.append(f'findings_path   = "{seed_csv}"')
            if findings_reference_median is not None:
                new_lines.append(
                    f'findings_reference_median = {findings_reference_median}'
                    '   # shared denominator (both samples) → preserves A/B ratios (PR #407)')
            if intensity_multiplier is not None:
                new_lines.append(
                    f'intensity_multiplier      = {intensity_multiplier}'
                    '   # brightness calibration; SAME for A and B under shared reference median')
        elif key == "save_path":
            new_lines.append(f'save_path       = "{save_path}"')
        elif key == "experiment_name":
            new_lines.append(f'experiment_name      = "{experiment_name}"')
        elif reference_path is not None and key == "reference_path":
            new_lines.append(f'reference_path  = "{reference_path}"')
        elif gradient_length is not None and key == "gradient_length":
            new_lines.append(f'gradient_length      = {gradient_length}'
                             '          # matched to the reference plasma .d span')
        elif superimpose and key == "add_real_data_noise":
            new_lines.append('add_real_data_noise           = false'
                             '   # stage 2: superimpose instead of sampling blank noise')
        elif superimpose and key == "superimpose_on_reference":
            new_lines.append('superimpose_on_reference      = true'
                             '    # overlay simulated spike-in on the real plasma frames')
        elif superimpose and key == "proteome_mix":
            # enforce the invariant: ratios stay baked in the seed, never re-diluted
            new_lines.append('proteome_mix = false')
        else:
            new_lines.append(line)
    return "\n".join(new_lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--seeds-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--samples", nargs="+", default=["A", "B"])
    ap.add_argument("--seed-suffix", default="",
                    help="Seed/variant suffix, e.g. '_imposed' to use "
                         "seed_sample<A|B>_imposed.csv and write to sample<A|B>_imposed/ "
                         "(keeps observed and imposed runs from colliding).")
    ap.add_argument("--path-mode", choices=["container", "host"], default="container")
    ap.add_argument("--repo-root", type=Path, default=None)
    ap.add_argument("--container-base", default="/work")
    ap.add_argument("--experiment-prefix", default="PLB-S1-DIA",
                    help="Experiment-name prefix (e.g. PLB-S2-DIA for the plasma stage).")
    ap.add_argument("--superimpose", action="store_true",
                    help="Stage 2: overlay simulated spike-in onto a real plasma .d "
                         "(add_real_data_noise=false, superimpose_on_reference=true). "
                         "Pair with --reference-path and YE-only seeds.")
    ap.add_argument("--reference-path", default=None,
                    help="Repoint reference_path (the real plasma .d for --superimpose).")
    ap.add_argument("--gradient-length", default=None,
                    help="Override gradient_length to match the reference plasma .d span.")
    ap.add_argument("--findings-reference-median", type=float, default=None,
                    help="Native ratio-preserving knob: shared event-scaling median for ALL "
                         "samples (sample A's median). Replaces the interim per-sample "
                         "--intensity-multiplier work-around (rustims PR #407).")
    ap.add_argument("--intensity-multiplier", type=float, default=None,
                    help="Brightness calibration written into the config; with "
                         "--findings-reference-median use ONE value for A and B.")
    args = ap.parse_args()

    repo_root = args.repo_root or Path(__file__).resolve().parents[1]
    base_text = args.base.read_text(encoding="utf-8")

    if args.path_mode == "host":
        subs = host_substitutions(repo_root)
        base_text = rewrite_paths(base_text, subs)
        work_prefix = subs["/work"]
    else:
        work_prefix = args.container_base

    rel_seed = args.seeds_dir.relative_to(repo_root) \
        if args.seeds_dir.is_absolute() else args.seeds_dir
    rel_save = args.out_dir.relative_to(repo_root) \
        if args.out_dir.is_absolute() else args.out_dir

    reference_path = args.reference_path
    if reference_path is not None and args.path_mode == "host":
        reference_path = rewrite_value(reference_path, host_substitutions(repo_root))

    suffix = args.seed_suffix
    for sample in args.samples:
        sample_dir = args.out_dir / f"sample{sample}{suffix}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        seed_path = f"{work_prefix}/{rel_seed}/seed_sample{sample}{suffix}.csv"
        save_path = f"{work_prefix}/{rel_save}/sample{sample}{suffix}"
        experiment = f"{args.experiment_prefix}-sample{sample}{suffix.replace('_', '-')}"
        text = render(base_text, sample=sample, seed_csv=seed_path,
                      save_path=save_path, experiment_name=experiment,
                      reference_path=reference_path, superimpose=args.superimpose,
                      gradient_length=args.gradient_length,
                      findings_reference_median=args.findings_reference_median,
                      intensity_multiplier=args.intensity_multiplier)
        out_path = sample_dir / "config.toml"
        out_path.write_text(text, encoding="utf-8")
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
