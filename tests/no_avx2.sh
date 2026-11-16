#!/bin/sh
# The run-time dispatch on CPUs without AVX2, under qemu-user CPU models: tests/no_avx2.sh [TREE]
# Builds yama's fuzz driver from TREE as the normal dispatch build and as a -DYAMA_SHADOW build, then on gen.py
# cases:
#   1. the dispatch build under `qemu-x86_64 -cpu Westmere` (no AVX2) must match a native -DYAMA_SIMD=0 build
#   2. the shadow build under Westmere must report p1_vec_rows=0 (the scalar path was taken), and under
#      `-cpu Haswell` (AVX2) p1_vec_rows>0, a positive control that the counter sees the AVX2 path
#   3. tests/run.sh with every program (the normal dispatch build) started under -cpu Westmere
# Needs qemu-x86_64 (qemu-user) and python3. Exits 1 on any FAIL.
set -u
H=$(cd "$(dirname "$0")" && pwd)
T=$(cd "${1:-$H/..}" && pwd) || exit 2
F=$H/yama_fuzz
command -v qemu-x86_64 > /dev/null || { echo "SKIP no_avx2 (no qemu-x86_64)"; exit 77; }
W=$(mktemp -d "${TMPDIR:-/tmp}/no-avx2.XXXXXX") || exit 2
trap 'rm -rf "$W"' EXIT
fail=0
b() { gcc -O2 -Wno-free-nonheap-object $2 -I"$T" "$F/main.c" "$T/mz_yama.c" "$T/mz_scores.c" "$T/util.c" -o "$W/$1"; }
b ref -DYAMA_SIMD=0 && b simd "" && b shadow -DYAMA_SHADOW || { echo "FAIL no_avx2: build"; exit 1; }
for s in 1 2 301; do
	g=""; [ $s = 301 ] && g="YAMA_FUZZ_GOP=268435456"
	python3 "$F/gen.py" $s 200 150 > "$W/in.$s"
	env $g "$W/ref" < "$W/in.$s" > "$W/ref.$s"
	env $g qemu-x86_64 -cpu Westmere "$W/simd" < "$W/in.$s" > "$W/qw.$s"
	if [ -s "$W/ref.$s" ] && cmp -s "$W/ref.$s" "$W/qw.$s"; then echo "PASS no_avx2 dispatch under Westmere = scalar, seed $s"
	else echo "FAIL no_avx2 dispatch under Westmere, seed $s"; fail=1; fi
done
# shadow counters: which path ran
for cpu in Westmere Haswell; do
	rm -f "$W/log.$cpu"
	YAMA_SHADOW_LOG="$W/log.$cpu" qemu-x86_64 -cpu $cpu "$W/shadow" < "$W/in.1" > /dev/null 2> "$W/err.$cpu"
	v=$(grep -o 'p1_vec_rows=[0-9]*' "$W/log.$cpu" "$W/err.$cpu" 2> /dev/null | head -1 | cut -d= -f2)
	echo "INFO no_avx2 shadow under $cpu: p1_vec_rows=${v:-?}"
	grep -q MISMATCH "$W/log.$cpu" "$W/err.$cpu" 2> /dev/null && { echo "FAIL no_avx2 shadow MISMATCH under $cpu"; fail=1; }
	if [ $cpu = Westmere ]; then [ "${v:-x}" = 0 ] || { echo "FAIL no_avx2: vector rows under Westmere"; fail=1; }
	else [ "${v:-0}" -gt 0 ] 2> /dev/null || { echo "FAIL no_avx2: positive control saw no vector rows under Haswell"; fail=1; }; fi
done
[ $fail = 0 ] && echo "PASS no_avx2 shadow: scalar path under Westmere, AVX2 path under Haswell"
# the whole regression suite (tests/run.sh) with every program started under qemu -cpu Westmere
make -s -C "$T" > /dev/null 2>&1 || { echo "FAIL no_avx2: make"; exit 1; }
mkdir -p "$W/bin"
for prog in roast multiz maf_project maf_sort maf_order tba multic all_bz pair2tb single_cov2 get_covered \
	get_standard_headers maf2fasta maf2lav lav2maf mafFind maf_checkThread blastzWrapper; do
	[ -x "$T/$prog" ] || continue
	printf '#!/bin/sh\nexec qemu-x86_64 -cpu Westmere "%s" "$@"\n' "$T/$prog" > "$W/bin/$prog"
	chmod +x "$W/bin/$prog"
done
if TMPDIR="$W" sh "$H/run.sh" "$W/bin" > "$W/runsh.out" 2>&1; then
	echo "PASS no_avx2 tests/run.sh under Westmere: $(tail -1 "$W/runsh.out")"
else
	echo "FAIL no_avx2 tests/run.sh under Westmere: $(tail -1 "$W/runsh.out")"; grep FAIL "$W/runsh.out" | head; fail=1
fi
exit $fail
