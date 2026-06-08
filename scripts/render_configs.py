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


def render(base_toml: str, *, sample: str, seed_csv: str, save_path: str,
           experiment_name: str) -> str:
    out = base_toml.replace("SAMPLE_PLACEHOLDER", sample)
    new_lines = []
    for line in out.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("findings_path"):
            new_lines.append(f'findings_path   = "{seed_csv}"')
        elif stripped.startswith("save_path"):
            new_lines.append(f'save_path       = "{save_path}"')
        elif stripped.startswith("experiment_name"):
            new_lines.append(f'experiment_name      = "{experiment_name}"')
        else:
            new_lines.append(line)
    return "\n".join(new_lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--seeds-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--samples", nargs="+", default=["A", "B"])
    ap.add_argument("--path-mode", choices=["container", "host"], default="container")
    ap.add_argument("--repo-root", type=Path, default=None)
    ap.add_argument("--container-base", default="/work")
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

    for sample in args.samples:
        sample_dir = args.out_dir / f"sample{sample}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        seed_path = f"{work_prefix}/{rel_seed}/seed_sample{sample}.csv"
        save_path = f"{work_prefix}/{rel_save}/sample{sample}"
        experiment = f"PLB-S1-DIA-sample{sample}"
        text = render(base_text, sample=sample, seed_csv=seed_path,
                      save_path=save_path, experiment_name=experiment)
        out_path = sample_dir / "config.toml"
        out_path.write_text(text, encoding="utf-8")
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
