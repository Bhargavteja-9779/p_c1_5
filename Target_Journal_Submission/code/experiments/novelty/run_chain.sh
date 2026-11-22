#!/bin/bash
# Resumable chain of the remaining novelty runs (each script skips finished folds/stages). Re-launch after any restart:
#   nohup setsid experiments/novelty/run_chain.sh >/dev/null 2>&1 &
set -u
cd "$(dirname "$0")/../.."
LOG=../results/experiment_logs; ST=$LOG/novelty_status.txt
run() { grep -q "$1" $ST 2>/dev/null && return; python3 experiments/novelty/$2 config/full.yaml >> $LOG/${2%.py}_full.stdout 2>&1 && echo "$1 $(date -u +%FT%TZ)" >> $ST; }
run NORMAN_STAGE_DONE run_norman_xsc.py
run CASE_STAGE_DONE run_norman_case.py
run CAND_STAGE_DONE run_candidates.py
run PIVOT_STAGE_DONE run_xsc_pivots.py
echo "CHAIN_DONE $(date -u +%FT%TZ)" >> $ST
