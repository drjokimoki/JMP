#!/usr/bin/env bash
# Waits for the endogenous rerun (PID 89806) to finish, then starts the
# H=6 exogenous T-scaling experiment.
# T=80 is already saved; the checkpoint logic in the script will skip it.
# Usage: bash launch_h6_scaling.sh

set -euo pipefail

WAIT_PID=89806
OUTDIR="RSM_exogenous_results_H6_scaling"
SCRIPT="rsm_common_loading_grid_weights_general_distance_reviewfix_final.R"
LOG="$OUTDIR/run_log_h6_rerun.txt"

cd "$(dirname "$0")"

if kill -0 "$WAIT_PID" 2>/dev/null; then
  echo "[$(date)] Waiting for PID $WAIT_PID (endogenous rerun) to finish..."
  while kill -0 "$WAIT_PID" 2>/dev/null; do
    sleep 60
  done
  echo "[$(date)] PID $WAIT_PID finished."
else
  echo "[$(date)] PID $WAIT_PID is not running — starting immediately."
fi

echo "[$(date)] Launching H=6 T-scaling..."

RSM_RUN_PROFILE=paper \
RSM_RUN_EXPERIMENTS=scaling_T_m \
RSM_SCALING_H=6 \
RSM_OUTDIR="$OUTDIR" \
  Rscript "$SCRIPT" > "$LOG" 2>&1

echo "[$(date)] Done. Results in $OUTDIR/"
