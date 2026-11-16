#!/usr/bin/env bash
# Are the Ix / Cx / Cz FITS32 range checks redundant under REAL-SHAPED gap tables?
#   tests/simd_mutation/real_table_checks.sh WORKDIR
# multiz's gap table has every nonzero gap-open entry equal (mz_scores.c), so the mutation report argued that
# the Ix, Cx and Cz checks are implied by the Iy, Cy, Dx and Dz checks that remain. This tests that claim
# directly: with YAMA_FUZZ_GOP=n every nonzero entry becomes n (the real shape, any magnitude), and the
# fits_i_removed and fits_cd_removed mutants must equal the scalar reference for every n tried, from the
# shipped value up past the 32-bit threshold, on random, targeted and real-call inputs.
# Prints SAME or DIFF per (mutant, n); exits 1 on any DIFF. (Single-entry tables, YAMA_FUZZ_GOP_IDX, are not
# real-shaped; check.sh uses them precisely because they are the only way to make one check matter alone.)
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$H/../.." && pwd)
W=${1:?usage: real_table_checks.sh WORKDIR}
mkdir -p "$W"
FZ=$SRC/tests/yama_fuzz
cc() { gcc -O2 -Wno-free-nonheap-object $2 -I"$1" "$FZ/main.c" "$1/mz_yama.c" "$1/mz_scores.c" "$1/util.c" -o "$3"; }
cc "$SRC" -DYAMA_SIMD=0 "$W/ref"
for m in fits_i_removed fits_cd_removed; do
	d=$W/$m; rm -rf "$d"; mkdir -p "$d"; cp "$SRC"/*.c "$SRC"/*.h "$d/"
	python3 "$H/mutants.py" apply "$m" "$d/mz_yama.c" > /dev/null || { echo "NOT APPLIED $m"; exit 2; }
	cc "$d" "" "$W/$m.bin"
done
inputs() {
	for s in 1 2 3; do python3 "$FZ/gen.py" $s 300 200; done
	for mode in narrow ties edges; do python3 "$FZ/gen.py" 7 300 200 $mode; done
	cat "$FZ/real_cases.txt"
}
inputs > "$W/in"
rc=0
# 0 = the shipped table (no override); then 2^k for k = 12 .. 31 (past INT_MAX/2 per entry)
for g in 0 4096 65536 1048576 16777216 67108864 134217728 268435456 536870912 1073741824 2147483647; do
	envs=(); [ $g != 0 ] && envs=(YAMA_FUZZ_GOP=$g)
	env "${envs[@]}" "$W/ref" < "$W/in" > "$W/ref.out" 2>&1
	for m in fits_i_removed fits_cd_removed; do
		env "${envs[@]}" "$W/$m.bin" < "$W/in" > "$W/m.out" 2>&1
		if cmp -s "$W/ref.out" "$W/m.out"; then echo "SAME $m gop=$g"; else echo "DIFF $m gop=$g"; rc=1; fi
	done
done
exit $rc
