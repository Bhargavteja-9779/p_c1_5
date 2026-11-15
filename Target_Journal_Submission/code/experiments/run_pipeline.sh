#!/bin/bash
# Full experimental pipeline, sequential and resumable (per-fold checkpoints; cached models). Re-run after interruption.
set -u
cd "$(dirname "$0")/.."
LOG=../results/experiment_logs
mkdir -p $LOG
python3 experiments/run_main.py config/full.yaml >> $LOG/run_main_full.stdout 2>&1
python3 experiments/run_followups.py config/full.yaml robust smrt hilic sens >> $LOG/run_followups_full.stdout 2>&1
echo "PIPELINE_DONE $(date -u +%FT%TZ)" >> $LOG/pipeline_status.txt
