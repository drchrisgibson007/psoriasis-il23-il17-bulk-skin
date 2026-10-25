#!/usr/bin/env bash
# End-to-end rerun from the GEO inputs in data/ and rma/ (about 1 minute). Needs Python 3 with numpy, pandas, scipy, matplotlib.
# Re-normalization from CEL files (3.6 GB, R/Bioconductor) is optional: see data/raw/README.md.
set -e
cd "$(dirname "$0")"
rm -rf out out_rma out_rma_qc; mkdir -p out out_rma out_rma_qc outputs
python3 scripts/overlap_check.py > /dev/null      # cohort overlap from GEO metadata -> outputs/overlap_check.txt
(cd scripts && python3 analyze.py && python3 sens.py && python3 analyze_rma.py) > outputs/run_log.txt 2>&1
# QC-dropped sensitivity run: writes lesion scores and summary, then stops at the GSE117468 response section
# because three QC-dropped baseline arrays are absent (expected; that section is not part of this sensitivity analysis)
(cd scripts && QCDROP=1 python3 analyze_rma.py) > outputs/run_log_qc.txt 2>&1 || echo "note: QC-dropped run stopped after writing out_rma_qc/ (expected, see comment)"
(cd scripts && python3 sens_rma.py) >> outputs/run_log.txt 2>&1
same=0; diff=0
for pair in "out:out_deposited" "out_rma:out_rma" "out_rma_qc:out_rma_qc"; do
  new=${pair%%:*}; ref=outputs_reference/${pair##*:}
  while IFS= read -r f; do r=${f#$ref/}
    if cmp -s "$f" "$new/$r"; then same=$((same+1)); else diff=$((diff+1)); echo "DIFFER $new/$r"; fi
  done < <(find "$ref" -type f ! -name run_log.txt | sort)
done
echo "MATCH $same files; DIFFER $diff"
cmp -s outputs/overlap_check.txt outputs_reference/overlap_check.txt && echo "MATCH  overlap_check.txt" || echo "DIFFER overlap_check.txt"
