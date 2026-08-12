# Legacy: the original Python Freddie

The first, Python version of the single-cell pipeline, kept for reference. It's
been replaced by the C++ [`freddie-core`](../freddie-core) (BAM → split → tint →
canonize → ILP → GTF) and the [`sclr-putagene`](../sclr-putagene) workflow, and
it isn't built or wired into the current Snakemake workflows.

The paths in `Snakefile` / `config.yaml` are relative to the old layout, so run
these from inside `legacy/` if you ever need them.

| file | what it was |
|------|-------------|
| `Snakefile`, `config.yaml` | the end-to-end workflow |
| `freddie.py`               | the CLI |
| `freddie/`                 | the package: `split`, `segment`, `ilp`, `isoforms`, `annotate`, `plot` — now in `freddie-core` |
| `Snakemake-envs/`          | per-rule conda envs |
