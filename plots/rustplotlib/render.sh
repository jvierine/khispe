#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export HDF5_DIR="${HDF5_DIR:-/opt/anaconda3}"
export DYLD_FALLBACK_LIBRARY_PATH="$HDF5_DIR/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
cargo build --release --manifest-path plots/rustplotlib/Cargo.toml
source_h5="${1:-output/ri-sweep-short-2h-20261013/animations/ri008-radar-full.h5}"
render_dir="${2:-output/ri-sweep-short-2h-20261013/animations/rustplotlib}"
plots/rustplotlib/target/release/syisr-rustplotlib "$source_h5" "$render_dir"
for field in ne ion_up; do
 ffmpeg -v error -y -framerate 30 -i "$render_dir/$field/%04d.png" -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -c:v libx264 -threads 2 -crf 19 -pix_fmt yuv420p -movflags +faststart "$render_dir/$field.mp4"
done
