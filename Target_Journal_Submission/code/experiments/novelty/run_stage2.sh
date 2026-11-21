#!/bin/bash
# Stage 2 of the novelty revision: waits for XSC, then candidate lists (PubChemLite) and NORMAN external XSC.
set -u
cd "$(dirname "$0")/../.."
LOG=../results/experiment_logs
until grep -q XSC_STAGE_DONE $LOG/novelty_status.txt 2>/dev/null; do sleep 60; done
grep -q CANDIDATES_DONE $LOG/run_candidates_full.stdout 2>/dev/null || python3 experiments/novelty/run_candidates.py config/full.yaml >> $LOG/run_candidates_full.stdout 2>&1
echo "CAND_STAGE_DONE $(date -u +%FT%TZ)" >> $LOG/novelty_status.txt
if [ -f experiments/novelty/run_norman_xsc.py ]; then
  grep -q NORMAN_XSC_DONE $LOG/run_norman_xsc_full.stdout 2>/dev/null || python3 experiments/novelty/run_norman_xsc.py config/full.yaml >> $LOG/run_norman_xsc_full.stdout 2>&1
  echo "NORMAN_STAGE_DONE $(date -u +%FT%TZ)" >> $LOG/novelty_status.txt
fi
