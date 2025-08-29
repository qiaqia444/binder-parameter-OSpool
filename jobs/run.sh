#!/usr/bin/env bash
set -euo pipefail
mkdir -p output logs
exec julia --project=. run_standard.jl "$@"
