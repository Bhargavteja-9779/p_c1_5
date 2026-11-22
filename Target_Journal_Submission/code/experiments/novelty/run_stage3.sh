#!/bin/bash
# Stage 3: real-world NORMAN suspect-screening case study (after NORMAN XSC multipliers exist).
set -u
cd "$(dirname "$0")/../.."
LOG=../results/experiment_logs
until grep -q NORMAN_STAGE_DONE $LOG/novelty_status.txt 2>/dev/null; do sleep 60; done
grep -q NORMAN_CASE_DONE $LOG/run_norman_case_full.stdout 2>/dev/null || python3 experiments/novelty/run_norman_case.py config/full.yaml >> $LOG/run_norman_case_full.stdout 2>&1
echo "CASE_STAGE_DONE $(date -u +%FT%TZ)" >> $LOG/novelty_status.txt
