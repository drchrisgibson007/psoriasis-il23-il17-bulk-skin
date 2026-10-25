#!/bin/bash
# Download the full (unslimmed) GEO inputs used by analyze.py / analyze_rma.py into ../data_full/
# Then point load.py at it:  sed -i "s#'..','data'#'..','data_full'#" ../scripts/load.py
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); OUT="$HERE/../data_full"; mkdir -p "$OUT"; cd "$OUT"
for g in GSE13355 GSE30999 GSE14905 GSE117468 GSE121212 GSE54456; do
  n=${g:3}; s=GSE$((10#$n/1000))nnn
  curl -sfL --retry 3 -o ${g}_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/$s/$g/matrix/${g}_series_matrix.txt.gz"
done
curl -sfL --retry 3 -o GSE121212_readcount.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121212/suppl/GSE121212_readcount.txt.gz"
curl -sfL --retry 3 -o GSE54456_RPKM_samples.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54456/suppl/GSE54456_RPKM_samples.txt.gz"
ls -l
