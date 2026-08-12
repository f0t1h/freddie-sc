#!/usr/bin/env bash
# Fetch the reference data the pipeline needs into data/ (kept out of git):
#
#   data/refs/homo_sapiens.dna.fa            GRCh38 primary assembly (soft-masked)
#   data/refs/homo_sapiens.cdna.fa           GRCh38 cDNA          (simulation only)
#   data/refs/homo_sapiens.annot.gtf         GRCh38 gene annotation
#   data/refs/homo_sapiens.chr21.*           chr21-only subset (fast sim tests)
#   data/whitelists/3M-february-2018.txt     10x 3' v3 barcode whitelist
#   data/whitelists/737K-august-2016.txt     10x 3' v2 / 5' v2 barcode whitelist
#
# data/ and raw-data/ are git-ignored — this script (re)creates their contents.
# cDNA is only used by the simulation (tksm2 model building); analysis-only runs
# don't need it.
#
# Usage:
#   bash scripts/fetch_data.sh              # download everything (needs network)
#   bash scripts/fetch_data.sh refs         # only the Ensembl references
#   bash scripts/fetch_data.sh whitelists   # only the 10x whitelists
#   ENSEMBL_RELEASE=110 bash scripts/fetch_data.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFS="$ROOT/data/refs"
WL="$ROOT/data/whitelists"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-110}"
EBASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens"
EGTF="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens"
# 10x whitelists ship inside Cell Ranger with no upstream download, so they're
# mirrored here. Override with WLBASE=<url or dir base> to use your own copy.
WLBASE="${WLBASE:-https://raw.githubusercontent.com/f0t1h/3M-february-2018/master}"

what="${1:-all}"

dl() {  # url  dest   (skips if dest exists; decompresses .gz)
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then echo "  [skip] $dest exists"; return; fi
  echo "  [get ] $url"
  if [[ "$url" == *.gz ]]; then
    curl -fL --retry 3 "$url" | gzip -dc > "$dest"
  else
    curl -fL --retry 3 "$url" -o "$dest"
  fi
}

fetch_refs() {
  mkdir -p "$REFS"
  echo "== Ensembl GRCh38 (release $ENSEMBL_RELEASE) -> data/refs =="
  dl "$EBASE/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" "$REFS/homo_sapiens.dna.fa"
  dl "$EBASE/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"               "$REFS/homo_sapiens.cdna.fa"
  dl "$EGTF/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"          "$REFS/homo_sapiens.annot.gtf"

  echo "== chr21 subset (fast sim tests) =="
  if command -v samtools >/dev/null 2>&1; then
    for kind in dna cdna; do
      if [[ ! -s "$REFS/homo_sapiens.chr21.$kind.fa" ]]; then
        samtools faidx "$REFS/homo_sapiens.$kind.fa" 21 \
          > "$REFS/homo_sapiens.chr21.$kind.fa" 2>/dev/null \
          || echo "  [warn] no '21' contig in homo_sapiens.$kind.fa (cdna uses tx ids; skipping)"
      fi
    done
  else
    echo "  [warn] samtools not found; skipping chr21 FASTA subset"
  fi
  # GTF chr21 subset: keep header comments + lines whose seqname is 21.
  if [[ ! -s "$REFS/homo_sapiens.chr21.annot.gtf" ]]; then
    awk -F'\t' '/^#/ || $1=="21"' "$REFS/homo_sapiens.annot.gtf" \
      > "$REFS/homo_sapiens.chr21.annot.gtf"
  fi
}

fetch_whitelists() {
  mkdir -p "$WL"
  echo "== 10x barcode whitelists -> data/whitelists =="
  dl "$WLBASE/3M-february-2018.txt.gz" "$WL/3M-february-2018.txt"   # 10x 3' v3
  dl "$WLBASE/737K-august-2016.txt"    "$WL/737K-august-2016.txt"   # 10x 3' v2 / 5' v2
}

case "$what" in
  refs)        fetch_refs ;;
  whitelists)  fetch_whitelists ;;
  all)         fetch_refs; fetch_whitelists ;;
  *) echo "usage: $0 [all|refs|whitelists]" >&2; exit 2 ;;
esac

echo "== done =="
