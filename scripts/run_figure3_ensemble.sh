#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p output/figure3-ensemble

pids=()
for config in configs/figure3-ensemble/*.json; do
    name="$(basename "$config" .json)"
    out="output/figure3-ensemble/$name"
    mkdir -p "$out"
    (
        unset DISPLAY
        mpirun -n 4 target/release/syisr-khi "$config" "$out" \
            >"$out/run.log" 2>&1
    ) &
    pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        failed=1
    fi
done
exit "$failed"
