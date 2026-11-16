#!/usr/bin/env bash
# Run multiz under a non-glibc allocator, to catch memory errors glibc's layout hides.
#   tests/allocator_check.sh BINDIR [N] [LEFT.maf RIGHT.maf]
# Runs BINDIR/multiz M=1 on a deep merge (default: bench/data/wfi_profile/in left.maf and right.maf) once under glibc,
# then N times (default 10) under jemalloc with plain settings and N times with MALLOC_CONF=junk:true (fresh and
# freed memory filled with junk). It reports each mode's crashes (exit status ≥ 128, e.g. 139 = SIGSEGV) and how many
# completed runs differ from the glibc output. Digests are of out.maf, '#' lines dropped.
# Exits 0 only if no run crashed and every run matched glibc; 77 if jemalloc is not installed.
#
# Why: pre_yama read map4[-1] and one byte past a column (mz_preyama.c lines 325 and 284, in multiz since the
# 2019 import). glibc puts a zeroed size-field half before every small chunk, so the stray read gave 0 and the
# output never showed it. jemalloc (like macOS malloc, which bioconda's osx builds use) has no such header, and
# stock 11.9 segfaulted in 5 of 33 deep merges.
set -uo pipefail
B=$(cd "${1:?usage: allocator_check.sh BINDIR [N] [LEFT RIGHT]}" && pwd)
N=${2:-10}
DEEP=/home/rcb112/project/claude/bench/data/wfi_profile/in
L=${3:-$DEEP/left.maf} R=${4:-$DEEP/right.maf}
JE=${JEMALLOC:-}
for c in /usr/lib/x86_64-linux-gnu/libjemalloc.so.2 /usr/lib/aarch64-linux-gnu/libjemalloc.so.2 /usr/local/lib/libjemalloc.so.2; do
	[ -z "$JE" ] && [ -e "$c" ] && JE=$c
done
[ -n "$JE" ] || { echo "SKIP allocator check: libjemalloc.so.2 not found (apt install libjemalloc2, or set JEMALLOC)"; exit 77; }
W=$(mktemp -d "${TMPDIR:-/tmp}/alloc-check.XXXXXX") || exit 2
trap 'rm -rf "$W"' EXIT
run() { # ENV... : one deep merge; prints "rc digest"
	local rc=0
	env "$@" "$B/multiz" M=1 "$L" "$R" 0 "$W/U1.maf" "$W/U2.maf" > "$W/out.maf" 2> /dev/null || rc=$?
	echo "$rc $(grep -v '^#' "$W/out.maf" 2> /dev/null | sha256sum | cut -c1-16)"
}
read -r grc ref < <(run)
echo "glibc: rc=$grc digest $ref"
[ "$grc" = 0 ] || { echo "FAIL: the glibc run itself failed"; exit 1; }
bad=0
for mode in plain junk; do
	envs=(LD_PRELOAD="$JE"); [ $mode = junk ] && envs+=(MALLOC_CONF=junk:true)
	crash=0 diff=0 ok=0 codes=""
	for ((i = 1; i <= N; i++)); do
		read -r rc dg < <(run "${envs[@]}")
		if [ "$rc" -ge 128 ] || [ "$rc" != 0 ]; then crash=$((crash + 1)); codes="$codes $rc"
		elif [ "$dg" != "$ref" ]; then diff=$((diff + 1))
		else ok=$((ok + 1)); fi
	done
	echo "jemalloc $mode: $N runs, $crash crashed${codes:+ (rc$codes)}, $diff differ from glibc, $ok identical"
	[ $crash = 0 ] && [ $diff = 0 ] || bad=1
done
[ $bad = 0 ] && echo "PASS allocator check" || echo "FAIL allocator check"
exit $bad
