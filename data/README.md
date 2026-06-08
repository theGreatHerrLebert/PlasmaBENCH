# data/

Raw + intermediate inputs. **Gitignored** — drop files at these repo-relative
paths (or point `PLB_DATA` at external storage; see `docker/.env.example`).

```
data/raw/
├── dia/blanks/   – TimsTOF DIA-PASEF blank(s) for noise injection (from Ute)
├── dia/plasma/   – real plasma-only .d (stage 2 superimpose target)
└── fasta/        – plasma + yeast + E. coli reference FASTA(s)
```

Inside the container these resolve under `/data/raw/...`. Configs reference
`/data/...`, never host paths, so they stay portable across machines.
