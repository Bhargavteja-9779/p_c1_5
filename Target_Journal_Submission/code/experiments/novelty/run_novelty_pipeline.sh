#!/bin/bash
# Novelty revision pipeline (Addendum G). Resumable: each stage skips finished folds/stages.
set -u
cd "$(dirname "$0")/../.."
LOG=../results/experiment_logs; mkdir -p $LOG
until grep -q PREWARM_DONE $LOG/prewarm.stdout 2>/dev/null; do
  pgrep -f novelty/prewarm.py >/dev/null || python3 experiments/novelty/prewarm.py config/full.yaml >> $LOG/prewarm.stdout 2>&1
  sleep 30
done
python3 experiments/novelty/run_xsc.py config/full.yaml >> $LOG/run_xsc_full.stdout 2>&1
echo "XSC_STAGE_DONE $(date -u +%FT%TZ)" >> $LOG/novelty_status.txt
