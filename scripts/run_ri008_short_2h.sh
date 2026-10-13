#!/usr/bin/env bash
set -Eeuo pipefail
cd /home/j/src/syisr-khi-rust
root=output/ri-sweep-short-2h-20261013
mkdir -p "$root/ri008"
printf '{"state":"queued","case":"ri008","waiting_for":"khi2d-ri-sweep-short-2h-20261013.service","mpi_ranks":128}\n' > "$root/ri008/status.json"
while systemctl --user is-active --quiet khi2d-ri-sweep-short-2h-20261013.service; do sleep 5; done
grep -q '"state":"complete"' "$root/status.json"
cfg=configs/ri-sweep-short-2h/ri008.json
out="$root/ri008"
test ! -f "$out/frame_00000.h5"
cp "$cfg" "$out/requested_config.json"
printf '{"state":"running","active_case":"ri008","stage":"simulation","mpi_ranks":128,"end_s":7200}\n' > "$root/status.json.tmp"
mv "$root/status.json.tmp" "$root/status.json"
trap 'printf "{\"state\":\"failed\",\"active_case\":\"ri008\"}\n" > "$root/status.json"' ERR
/usr/bin/mpirun -np 128 --bind-to core target/release/syisr-khi "$cfg" "$out" > "$out/run.log" 2>&1
test -s "$out/frame_00720.h5"
test -s "$out/diagnostics.json"
target/release/pack-density "$out" "$out/density_frames.h5" > "$out/packing.log" 2>&1
printf '{"case":"ri008","state":"complete","completed_utc":"%s"}\n' "$(date -u +%FT%TZ)" > "$out/complete.json"
printf '{"state":"complete","active_case":"none","stage":"finished","case_count":6,"mpi_ranks":128,"end_s":7200}\n' > "$root/status.json.tmp"
mv "$root/status.json.tmp" "$root/status.json"
cp "$out/complete.json" "$out/status.json"
