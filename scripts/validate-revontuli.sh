#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
export PATH="$HOME/.cargo/bin:$PATH"
export OMP_NUM_THREADS=1
unset DISPLAY XAUTHORITY
mkdir -p output/validation
printf '%s\n' '{"nx":32,"nz":48,"validate_only":true}' > output/validation/coarse.json
printf '%s\n' '{"nx":64,"nz":96,"validate_only":true}' > output/validation/fine.json
for ranks in 1 2 4; do
  mpirun -np "$ranks" target/release/syisr-khi output/validation/coarse.json "output/validation/poisson-$ranks"
  mpirun -np "$ranks" target/release/syisr-khi configs/smoke.json "output/validation/ranks-$ranks"
done
mpirun -np 4 target/release/syisr-khi output/validation/fine.json output/validation/poisson-fine
mpirun -np 2 target/release/syisr-khi configs/neutral-control.json output/validation/control
target/release/syisr-khi --compare output/validation/ranks-1/frame_00001.h5 output/validation/ranks-2/frame_00001.h5
target/release/syisr-khi --compare output/validation/ranks-1/frame_00001.h5 output/validation/ranks-4/frame_00001.h5
target/release/check-validation
