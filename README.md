# scFreddie

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg?style=flat)](https://snakemake.github.io)

- [scFreddie](#scfreddie)
	- [Running scFreddie](#running-scfreddie)
		- [Using Snakemake](#using-snakemake)
		- [Running the tools directly](#running-the-tools-directly)
		- [Containers](#containers)
	- [Configuration](#configuration)
		- [Samplesheet](#samplesheet)
		- [Tool binaries](#tool-binaries)
		- [Apptainer profile](#apptainer-profile)
	- [External data](#external-data)
	- [Pipeline stages](#pipeline-stages)
	- [Simulation](#simulation)
	- [Solver](#solver)
- [Citing scFreddie](#citing-scfreddie)

scFreddie extends [Freddie](https://github.com/vpc-ccg/freddie)—an annotation-independent tool for identifying and discovering isoforms from long-read transcriptomic data—to single-cell long reads. It calls isoforms per cell without requiring the splice aligner (e.g. `minimap2`) to be supplied with known isoform annotations, taking raw single-cell long reads through to a per-cell isoform `GTF`.

The repository holds two independent Snakemake workflows:

- **Analysis** ([`sclr-putagene/`](sclr-putagene)): the pipeline itself—barcode demultiplexing, alignment, per-cell locus construction, cell clustering, isoform calling, and evaluation.
- **Simulation** (`Simulate.smk` + [`extern/tksm2/`](extern/tksm2)): generates 10x single-cell long reads with known ground truth for benchmarking the analysis. Documented separately in [`Simulate.README.md`](Simulate.README.md).

The isoform caller invoked by the analysis lives in [`freddie-core/`](freddie-core) (a C++ rewrite; `BAM` → split → tint → canonize → ILP → `GTF`). The original Python Freddie is preserved under [`legacy/`](legacy) for reference. Config paths are all relative to the repository and external data is fetched by a script, so nothing here is pinned to one machine.

## Running scFreddie

The two workflows are separate—there is no shared entry point. Run the analysis on its own; run the simulator first only when the samplesheet points at simulated FASTQs (`output/simulate/TS/*.fastq`), which the simulator produces.

### Using Snakemake

Clone with submodules, create the environment, fetch the reference data, and build the C++ tools:

```bash
git clone --recursive <repo-url> scfreddie && cd scfreddie
# forgot --recursive?  git submodule update --init --recursive

micromamba create -f environment.yaml && micromamba activate scfreddie
bash scripts/fetch_data.sh          # genome, annotation, barcodes -> data/

cmake -S . -B build                 # builds freddie-core, putagene, tksm2
cmake --build build -j              # (-DFREDDIE_ENABLE_GUROBI=ON for Gurobi)
```

Then run the analysis (and the simulator first, if your samplesheet references simulated FASTQs):

```bash
snakemake -s sclr-putagene/Snakefile --configfile config.yaml         -j8
```

### Running the tools directly

The compiled tools are ordinary command-line programs and can be run outside Snakemake:

```bash
putagene --help       # per-cell locus construction and read grouping
freddie-core --help   # BAM -> isoform GTF
```

### Containers

The container image is a **tool bundle**: it compiles `freddie-core` and `putagene` on top of a base that already carries the conda dependencies.

```bash
# With Docker
docker build -f docker/Dockerfile -t scfreddie:latest . 
docker run --rm scfreddie:latest putagene --help
docker run --rm scfreddie:latest freddie-core --help

# With apptainer (singularity)
apptainer build scfreddie.sif docker/scfreddie.def 
apptainer exec scfreddie.sif putagene --help
apptainer exec scfreddie.sif freddie-core --help
```

## Configuration

### Samplesheet

Samples live in a CSV, `samplesheet.csv`, which keeps the sample list out of the tool config (`config.yaml`) and the simulator config (`Simulate.config.yaml`). One row per sample:

```csv
sample,fastq,ref,chemistry,whitelist,demuxed
NOA_chr19,raw-data/NOA_chr19.fastq.gz,GRCh38,10x3v3,data/whitelists/3M-february-2018.txt,true
S1_3v3_recal,output/simulate/TS/S1_3v3_recal.fastq,GRCh38,10x3v3,data/whitelists/3M-february-2018.txt,false

```
Required columns:
- `sample`: sample name; used for output paths.
- `fastq`: a real file under `raw-data/` or a simulated one under `output/simulate/TS/` (Snakemake builds the simulated ones on demand from `Simulate.config.yaml`).
- `ref`: names an entry in the config's `refs`.
- `chemistry`: 10x chemistry, e.g. `10x3v3`.
- `whitelist`: 10x barcode whitelist path.
- `demuxed`: `true` when the reads already carry barcodes and should skip flexiplex.

Point at a different sheet with `--config samplesheet=other.csv`.

### Tool binaries

The analysis resolves `freddie-core` and `putagene` in order: the config's `binaries:` block (an absolute path or a bare command), then `PATH`, otherwise it stops with a clear error. After `cmake --build build`, either put `build/` on `PATH` or point the config at the built binaries:

```yaml
# config.yaml
binaries:
    freddie-core: /abs/path/freddie-core/build/freddie-core
    putagene:     /abs/path/sclr-putagene/build/putagene
```

The workflow does not build the tools itself, and it does not assume it is run from the repository root—its own scripts are found relative to the Snakefile, so `snakemake -s /abs/sclr-putagene/Snakefile ...` works from anywhere. Data paths (`outpath`, `samplesheet`, `refs`) are yours to set; make them absolute when running outside the repository.

### Apptainer profile

`profiles/apptainer/` runs every rule inside `scfreddie.sif` (the tools are on the image's `PATH`). Build the image, edit the profile's `.sif` path and `--bind`, then:

```bash
apptainer build scfreddie.sif docker/scfreddie.def
snakemake --workflow-profile profiles/apptainer \
          -s sclr-putagene/Snakefile --configfile config.yaml -j8
```

## Fetching genome reference and barcodes

`scripts/fetch_data.sh` fills the (git-ignored) `data/` tree:

```
data/refs/homo_sapiens.{dna,cdna}.fa, .annot.gtf         Ensembl GRCh38
data/refs/homo_sapiens.chr21.*                           chr21 only, for quick runs
data/whitelists/{3M-february-2018,737K-august-2016}.txt  10x barcodes
```

The genome, cDNA, and annotation come from Ensembl (cDNA is only needed for
simulation). The 10x whitelists ship inside Cell Ranger with no upstream download,
so they're pulled from a [mirror](https://github.com/f0t1h/3M-february-2018);
point at your own copy with `WLBASE=<url-or-dir-base> bash scripts/fetch_data.sh
whitelists`.

Real reads go in `raw-data/`, pipeline output in `output/`; both are git-ignored.

## Pipeline stages

The analysis accepts single-cell long reads and outputs a per-cell isoform `GTF`, wired together by Snakemake through these stages:

- **flexiplex**: detects and extracts cell barcodes and UMIs from the raw reads (skipped when `demuxed` is `true`).
- **minimap2**: splice-aware alignment of the demultiplexed reads to the reference genome.
- **putagene**: per-cell locus construction and read grouping (deduplication, strand handling).
- **Scanpy**: cell clustering from the count matrix, with optional Scrublet doublet removal.
- **freddie-core**: annotation-independent isoform calling per cluster (`BAM` → split → tint → canonize → ILP → `GTF`).
- **eval**: scores the called isoforms against the ground truth (simulated samples only).

## Simulation

Its configuration, pipeline modules, and outputs are documented in [`Simulate.README.md`](Simulate.README.md).

## Solver

freddie-core solves with HiGHS by default, and it is what the containers use. Gurobi is optional for **local** builds on a machine with gurobi installed and configure with `-DFREDDIE_ENABLE_GUROBI=ON`.

# Citing scFreddie

scFreddie builds on Freddie, which was accepted to [RECOMB 2021](http://web.archive.org/web/20220129112349/https://www.recomb2021.org/program) and published in Nucleic Acids Research. If you use scFreddie in your research, please cite:

*Freddie: Annotation-independent Detection and Discovery of Transcriptomic Alternative Splicing Isoforms*. Baraa Orabi, Brian McConeghy, Ning Xie, Xuesen Dong, Cedric Chauve, Faraz Hach. Nucleic Acids Research 2022; DOI: [10.1093/nar/gkac1112](https://doi.org/10.1093/nar/gkac1112)

