#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
export PATH="$HOME/.cargo/bin:$PATH"
export MPICC="$PWD/scripts/mpicc-local"
cargo test --locked
cargo build --release --locked
