# Simulation (`Simulate.smk`)

Generates synthetic single-cell long-read data with known ground truth, so the
analysis pipeline can be benchmarked against a real answer. It's a thin wrapper
over [tksm2](extern/tksm2) (`extern/tksm2/Snakefile` is imported as a module):
each experiment is a pipeline of tksm modules that builds molecules, gives them
10x barcodes, truncates and amplifies them, and finally sequences them with
badread-style error. The output is a FASTQ plus a ground-truth TSV per experiment.

This workflow is independent of the analysis — run it on its own. Its FASTQs feed
the analysis only through the samplesheet (see below).

```sh
snakemake -s Simulate.smk --configfile Simulate.config.yaml -j8
```

Outputs land under `output/simulate/` (config `outpath`):

```
output/simulate/TS/<experiment>.fastq   simulated reads
output/simulate/TS/<experiment>.tsv     ground truth (per-read origin: cell, transcript, ...)
output/simulate/preprocess/             cached tksm intermediates
```

## Config

`Simulate.config.yaml` has five parts:

- **`TS_experiments`** — the experiments to build. Each is a named `pipeline:` — an
  ordered list of tksm modules, each with `params:` (raw tksm flags) and sometimes
  a `model:` (a named entry under `models:`). Modules chain: one's output MDF is
  the next one's input.
- **`refs`** — genome/annotation/cDNA sets (`GTF`/`DNA`/`cDNA`) and barcode
  whitelists, pulled from `data/` (see *External data* in the main [README](README.md)).
- **`samples`** — the real reference reads the expression models are learned from
  (SGNex direct-cDNA runs under `data/samples/`), each tied to a `ref`.
- **`models`** — per-experiment tksm model definitions (currently `Tsb`
  transcribe/sample models: which `sample`, which barcode list `cb-txt`, and flags
  like `--cb-count`, `--cb-lognorm-params`).
- **`exec`** — `tksm: extern/tksm2/build/tksm2`, the built simulator binary.

### Pipeline modules

The modules used across the experiments, in the order they typically appear:

| module | what it does |
|--------|--------------|
| `Tsb` | transcribe & sample molecules from a real expression profile (per cell barcode) |
| `Mrg` | merge several sources (the cell types C1–C5) into one population |
| `plA` | add a poly-A tail |
| `Flp` | strand-flip a fraction of molecules |
| `Tag` | prepend/append fixed sequence (adapters, primers) or a random N-mer (UMI) |
| `SCB` | attach the single-cell barcode |
| `Trc` | truncate reads (models real 3'/5' coverage loss) |
| `PCR` | amplify with a per-cycle error model |
| `Shf` | shuffle molecule order |
| `Seq` | sequence to FASTQ with badread-style error |

The `C1`–`C5` experiments are the five cell types (A549, Hct116, HepG2, K562,
MCF7). The `S1*` experiments merge those into a 10x library of a given chemistry
(`_3v3`, `_5v2`, `_3v2`); the barcode-bearing end is built before `Trc` so
truncation preserves it. Variants like `_recal` / `_alwaysend` / `_notrunc` only
change the truncation step — see the comments inline in the config for why.

## Feeding the analysis

The analysis consumes a simulated experiment by referencing its FASTQ in the
samplesheet (`demuxed: false`, so it goes through barcode demux like a real
sample):

```csv
sample,fastq,ref,chemistry,whitelist,demuxed
S1_3v3_recal,output/simulate/TS/S1_3v3_recal.fastq,GRCh38,10x3v3,data/whitelists/3M-february-2018.txt,false
```

Snakemake builds `output/simulate/TS/S1_3v3_recal.fastq` on demand from the
matching `TS_experiments` entry. The eval rule also reads the ground-truth
`output/simulate/TS/<experiment>.tsv` to score the analysis output.

## Data

Beyond the shared references and whitelists in `data/` (main README), simulation
needs the real reference reads the models learn from, under `data/samples/` — the
SGNex direct-cDNA FASTQs named in the config's `samples:` block. tksm2 must be
built (`extern/tksm2/build/tksm2`, from the top-level `cmake --build build`).
