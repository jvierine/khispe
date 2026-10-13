#!/usr/bin/env bash
set -Eeuo pipefail
cd /home/j/src/syisr-khi-rust
root=output/ri-sweep-short-2h-20261013
mkdir -p "$root"
cp configs/ri-sweep-short-2h/manifest.json "$root/manifest.json"
sha256sum target/release/syisr-khi > "$root/solver.sha256"
active=none
stage=start
write_status() {
    local state=$1
    printf '{"state":"%s","active_case":"%s","stage":"%s","updated_utc":"%s","mpi_ranks":128,"end_s":7200}\n' \
        "$state" "$active" "$stage" "$(date -u +%FT%TZ)" > "$root/status.json.tmp"
    mv "$root/status.json.tmp" "$root/status.json"
}
trap 'write_status failed' ERR
trap 'write_status stopped; exit 130' INT TERM
for active in ri010 ri012 ri014 ri016 ri018; do
    out="$root/$active"
    cfg="configs/ri-sweep-short-2h/$active.json"
    mkdir -p "$out"
    if [[ -f "$out/frame_00000.h5" ]]; then
        echo "Refusing to overwrite existing run: $out" >&2
        exit 1
    fi
    cp "$cfg" "$out/requested_config.json"
    stage=simulation
    write_status running
    /usr/bin/mpirun -np 128 --bind-to core target/release/syisr-khi "$cfg" "$out" > "$out/run.log" 2>&1
    test -s "$out/frame_00720.h5"
    test -s "$out/diagnostics.json"
    stage=packing
    write_status running
    target/release/pack-density "$out" "$out/density_frames.h5" > "$out/packing.log" 2>&1
    printf '{"case":"%s","state":"complete","completed_utc":"%s"}\n' \
        "$active" "$(date -u +%FT%TZ)" > "$out/complete.json"
done
active=none
stage=finished
write_status complete
