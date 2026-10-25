#!/bin/bash
# usage: bash fetch_raw.sh DEST   (reads raw_manifest.tsv in this folder; needs curl, tar, md5sum)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); DEST=${1:-raw_download}; mkdir -p "$DEST"; cd "$DEST"
tail -n +2 "$HERE/raw_manifest.tsv" | cut -f1 | sort -u | while read acc; do
  case $acc in
    GSE*) n=${acc:3}; if [ ${#n} -le 3 ]; then s=GSEnnn; else s=GSE$((10#$n/1000))nnn; fi
          if grep -qP "^$acc\traw/$acc/GSM|^$acc\traw/${acc}_RAW.tar" "$HERE/raw_manifest.tsv"; then
            [ -s raw/${acc}_RAW.tar ] || { mkdir -p raw; curl -sfL --retry 3 -o raw/${acc}_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/$s/$acc/suppl/${acc}_RAW.tar"; }
            if grep -qP "^$acc\traw/$acc/" "$HERE/raw_manifest.tsv"; then mkdir -p raw/$acc && tar -xf raw/${acc}_RAW.tar -C raw/$acc; fi
          fi ;;
  esac
done
# any listed file still missing: fetch from its URL (GEO series files, ArrayExpress files, NCBI counts, series matrices)
tail -n +2 "$HERE/raw_manifest.tsv" | while IFS=$'\t' read acc rel bytes md5 url orig; do
  [ -s "$rel" ] && continue
  case "$url" in http*" (inside)"|"see notes") continue;; esac
  mkdir -p "$(dirname "$rel")"; curl -sfL --retry 3 -o "$rel" "$url" || echo "download failed: $rel"
done
echo "checking MD5..."; bad=0
tail -n +2 "$HERE/raw_manifest.tsv" | while IFS=$'\t' read acc rel bytes md5 url orig; do
  if [ ! -s "$rel" ]; then echo "MISSING $rel"; continue; fi
  [ "$(md5sum "$rel" | cut -d' ' -f1)" = "$md5" ] || echo "MD5 DIFFERS $rel"
done
echo "done (lines above starting MISSING or MD5 DIFFERS need attention)"
