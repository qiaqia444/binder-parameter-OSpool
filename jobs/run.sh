#!/usr/bin/env bash
set -euo pipefail
REPO_DIR="$(pwd)"
SIF_PATH="${SIF_PATH:-$REPO_DIR/containers/container.sif}"
APPT="${APPTAINER:-apptainer}"
mkdir -p "$REPO_DIR/output"
"$APPT" exec --bind "$REPO_DIR":/workdir "$SIF_PATH" \
  julia --project=/workdir -t ${JULIA_NUM_THREADS:-1} /workdir/run.jl \
    ${PARAMS_JSON:+--params "$PARAMS_JSON"} \
    --outdir /workdir/output
