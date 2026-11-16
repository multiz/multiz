#!/usr/bin/env bash
# Interleaved timing of multiz builds, one job at a time.
#   simd_row_timing.sh deep  BINROOT OUTDIR ROUNDS BUILD...   (cannabis deep merge, M=1)
#   simd_row_timing.sh piece BINROOT OUTDIR ROUNDS BUILD...   (hg38 100-way piece 03)
# BINROOT/BUILD holds roast, multiz, maf_project.  Each round runs every build once, in a rotated
# order.  Prints one TSV line per run: kind build round wall_s user_s sys_s digest (non-'#' lines).
set -euo pipefail
kind=$1 BR=$2 OUT=$3 R=$4; shift 4; builds=("$@")
DEEP=${DEEP:-/home/rcb112/project/claude/bench/data/wfi_profile/in}
U=${U:-/home/rcb112/project/claude/bench/data/ucsc100way}
PIECE=${PIECE:-/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts/ucsc100way_run_piece.sh}
mkdir -p "$OUT"
nb=${#builds[@]}
for ((r = 1; r <= R; r++)); do
  for ((i = 0; i < nb; i++)); do
    b=${builds[$(( (i + r - 1) % nb ))]}
    d=$OUT/$kind.$b.$r
    rm -rf "$d"; mkdir -p "$d"
    if [ "$kind" = deep ]; then
      /usr/bin/time -f "%e %U %S" -o "$d/time" "$BR/$b/multiz" M=1 "$DEEP/left.maf" "$DEEP/right.maf" 0 \
        "$d/U1.maf" "$d/U2.maf" > "$d/out.maf" 2> "$d/err"
      read -r wall user sys < "$d/time"
      dig=$(cat "$d/out.maf" "$d/U1.maf" "$d/U2.maf" | grep -v '^#' | sha256sum | cut -c1-16)
      rm -f "$d"/*.maf
    else
      bash "$PIECE" "$BR/$b" "$U/run03/inputs-v3" "$U/tree.nh" "$d/run" > /dev/null
      read -r _ wall _ _ _ _ user _ sys < "$d/run/elapsed.txt"
      dig=$(grep -v '^#' "$d/run/out.maf" | sha256sum | cut -c1-64)
      rm -rf "$d/run/tmp" "$d/run/out.maf"
    fi
    printf '%s\t%s\t%d\t%s\t%s\t%s\t%s\n' "$kind" "$b" "$r" "$wall" "$user" "$sys" "$dig"
  done
done
