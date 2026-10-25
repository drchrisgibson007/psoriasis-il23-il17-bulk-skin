#!/usr/bin/env bash
# Run from anywhere:  bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
source config.sh
bash analysis/run_all.sh
echo "run_all: finished"
