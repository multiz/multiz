#!/bin/sh
# make check's yama fuzz test: tests/yama_fuzz/check.sh [TREE]
# Builds yama from TREE (default: this source tree) twice, as a scalar reference (-DYAMA_SIMD=0) and as the
# normal run-time-dispatch build, and requires byte-identical output on random cases (gen.py):
#   plain       seeds 1-3, 300 cases
#   fits32      seeds 301-302 with YAMA_FUZZ_GOP=2^28, so row coefficients pass 32 bits and the AVX2 rows
#               must fall back to scalar
#   fits32 idx  seeds 401-402 with YAMA_FUZZ_GOP_IDX=2, 0 and 8: one FITS32 check at a time (the shipped
#               gap table makes the checks imply each other; see main.c). These are the only cases that
#               catch a removed range check (tests/simd_mutation).
#   fits32 wA   seeds 601-602 with YAMA_FUZZ_SS=2^23: substitution scores scaled so a row's summed scores (wA)
#               pass 32 bits while the value-bound check still allows AVX2; the only cases that catch a removed
#               wA FITS32 check
#   override    one plain seed with YAMA_SIMD=0 in the environment, which must force the scalar path
#   targeted    gen.py's narrow, ties and edges modes, seeds 501-503 each: band widths 1-9 and 10-17 (every
#               vector tail), MININT neighbours, forced C/D/I ties including D/I ties, rows 0 and M, the
#               first-column waiver and all-dash edge columns
#   bound       the value-bound check (YAMA_VBOUND, YAMA_ZBOUND, YAMA_EBOUND in mz_yama.c), on a plain seed with
#               -DYAMA_SHADOW builds: with the shipped limits it never fires (bound_fallbacks=0) and AVX2 rows run;
#               with each limit lowered to 0 in turn it fires, and the output still equals the scalar reference;
#               with YAMA_VBOUND=0 it fires on every call and no row reaches either AVX2 kernel. Needs AVX2 natively
#               (the check sits on the AVX2 path), so it is skipped under YAMA_FUZZ_RUN, without AVX2, or with
#               YAMA_SIMD=0 in the environment (tests/validate.sh runs the suite that way too)
#   real        tests/yama_fuzz/real_cases.txt: yama calls captured from multiz on tests/data/synth/k2l2s1 with a
#               -DYAMA_DUMP_CASES build (pre_yama's refinement call: one row of A against two of B). They are the
#               only fuzz input found to expose a D/I tie decided in the vector kernel (the tie_di mutant)
# Every comparison requires both drivers to exit 0 and compares stdout only (a crash in both builds, or the same
# message on stderr, used to pass). On a CPU without AVX2, or with YAMA_SIMD=0 in the environment, both builds run
# scalar and the comparison is vacuous; it says so. Needs python3.
# Prints PASS/FAIL lines and exits 1 on any FAIL, 0 otherwise (77 = skipped, no python3).
set -u
H=$(cd "$(dirname "$0")" && pwd)
T=$(cd "${1:-$H/../..}" && pwd) || exit 2
command -v python3 > /dev/null 2>&1 || { echo "SKIP yama fuzz (no python3)"; exit 77; }
W=$(mktemp -d "${TMPDIR:-/tmp}/yama-fuzz.XXXXXX") || exit 2
trap 'rm -rf "$W"' EXIT
trap 'exit 2' HUP INT TERM
GOP=268435456
# YAMA_FUZZ_CC and YAMA_FUZZ_CFLAGS (default gcc, -O2) choose the compiler and flags for both builds, so the
# sanitizer and compiler-matrix runs (tests/validate.sh) cover this check too. make check passes $(CC) and
# $(CFLAGS), so "make CC=clang check" fuzzes clang's build, and a machine without gcc can run it.
FCC=${YAMA_FUZZ_CC:-gcc}
FCFLAGS=${YAMA_FUZZ_CFLAGS:--O2}
# YAMA_FUZZ_RUN prefixes every run of the built drivers, e.g. qemu-aarch64 for a cross build (tests/validate.sh cross)
FRUN=${YAMA_FUZZ_RUN:-}
build() {
	$FCC $FCFLAGS -Wno-free-nonheap-object $2 -I"$T" "$H/main.c" "$T/mz_yama.c" "$T/mz_scores.c" "$T/util.c" \
		-o "$W/$1" || { echo "FAIL yama fuzz: cannot build $1"; exit 1; }
}
build ref -DYAMA_SIMD=0
build simd ""
AVX2=0
if [ -z "$FRUN" ] && grep -q -w avx2 /proc/cpuinfo 2> /dev/null; then AVX2=1; else
	echo "NOTE yama fuzz: this CPU has no AVX2, so both builds take the scalar path"
fi
if [ "${YAMA_SIMD-}" = 0 ]; then
	echo "NOTE yama fuzz: YAMA_SIMD=0 is set, so both builds take the scalar path"
fi
fail=0
# same LABEL A B: both runs exited 0 (A.rc, B.rc), and A's stdout is non-empty and equals B's
same() {
	if [ "$(cat "$2.rc")" = 0 ] && [ "$(cat "$3.rc")" = 0 ] && [ -s "$2" ] && cmp -s "$2" "$3"; then
		echo "PASS $1"
	else
		echo "FAIL $1 (exit $(cat "$2.rc") and $(cat "$3.rc"))"; fail=1
	fi
}
# drive BIN IN OUT [ENV...]: run a driver; stdout to OUT, stderr to OUT.err, exit status to OUT.rc
drive() {
	bin=$1 in=$2 out=$3; shift 3
	env "$@" $FRUN "$bin" < "$in" > "$out" 2> "$out.err"
	echo $? > "$out.rc"
}
# cmpcase LABEL SEED NCASES [ENV...]: run both builds on gen.py SEED NCASES 200 $MODE under ENV, compare
MODE=random
cmpcase() {
	label=$1 seed=$2 nc=$3; shift 3
	python3 "$H/gen.py" "$seed" "$nc" 200 "$MODE" > "$W/in" || { echo "FAIL $label seed $seed: gen.py"; fail=1; return; }
	drive "$W/ref" "$W/in" "$W/ref.out" "$@"
	drive "$W/simd" "$W/in" "$W/simd.out" "$@"
	same "$label seed $seed" "$W/ref.out" "$W/simd.out"
}
for s in 1 2 3; do cmpcase "yama fuzz plain" $s 300; done
for s in 301 302; do cmpcase "yama fuzz fits32" $s 300 YAMA_FUZZ_GOP=$GOP; done
for ix in 2 0 8; do
	for s in 401 402; do cmpcase "yama fuzz fits32 idx $ix" $s 300 YAMA_FUZZ_GOP=$GOP YAMA_FUZZ_GOP_IDX=$ix; done
done
for s in 601 602; do cmpcase "yama fuzz fits32 wA" $s 300 YAMA_FUZZ_SS=8388608; done
for MODE in narrow ties edges; do
	for s in 501 502 503; do cmpcase "yama fuzz targeted $MODE" $s 300; done
done
MODE=random
# the value-bound check: shadow builds report bound_fallbacks (calls sent to scalar by the check) at exit
shadow_field() { tr ' ' '\n' < "$1" | sed -n "s/^$2=//p" | head -1; }
if [ "$AVX2" = 1 ] && [ "${YAMA_SIMD-}" != 0 ]; then
	python3 "$H/gen.py" 9 300 200 > "$W/in"
	drive "$W/ref" "$W/in" "$W/ref.out"
	for b in ctl:"" v:-DYAMA_VBOUND=0 z:-DYAMA_ZBOUND=0 e:-DYAMA_EBOUND=0; do
		n=${b%%:*}
		build "b$n" "-DYAMA_SHADOW ${b#*:}"
		drive "$W/b$n" "$W/in" "$W/b$n.out"
		fb=$(shadow_field "$W/b$n.out.err" bound_fallbacks) calls=$(shadow_field "$W/b$n.out.err" calls)
		p1=$(shadow_field "$W/b$n.out.err" p1_vec_rows) p2=$(shadow_field "$W/b$n.out.err" p2_vec_rows)
		case $n in
		ctl) ok=$([ "${fb:-x}" = 0 ] && [ "${p1:-0}" -gt 0 ] && [ "${p2:-0}" -gt 0 ] && echo 1) ;;
		v) ok=$([ "${fb:-0}" -gt 0 ] && [ "$fb" = "$calls" ] && [ "${p1:-x}" = 0 ] && [ "${p2:-x}" = 0 ] && echo 1) ;;
		*) ok=$([ "${fb:-0}" -gt 0 ] && echo 1) ;;
		esac
		if [ "$ok" = 1 ] && [ "$(cat "$W/ref.out.rc")" = 0 ] && [ "$(cat "$W/b$n.out.rc")" = 0 ] \
				&& [ -s "$W/ref.out" ] && cmp -s "$W/ref.out" "$W/b$n.out"; then
			echo "PASS yama fuzz bound $n (bound_fallbacks=$fb of $calls calls, p1_vec_rows=$p1)"
		else
			echo "FAIL yama fuzz bound $n (bound_fallbacks=${fb:-?} of ${calls:-?} calls, p1_vec_rows=${p1:-?})"; fail=1
		fi
	done
	# the run-time override, seen from inside: with YAMA_SIMD=0 the dispatch build sends no row to AVX2
	drive "$W/bctl" "$W/in" "$W/ov.out" YAMA_SIMD=0
	p1=$(shadow_field "$W/ov.out.err" p1_vec_rows) p2=$(shadow_field "$W/ov.out.err" p2_vec_rows)
	if [ "${p1:-x}" = 0 ] && [ "${p2:-x}" = 0 ] && [ "$(cat "$W/ov.out.rc")" = 0 ] && cmp -s "$W/ref.out" "$W/ov.out"; then
		echo "PASS yama fuzz override YAMA_SIMD=0 keeps every row scalar"
	else
		echo "FAIL yama fuzz override YAMA_SIMD=0 (p1_vec_rows=${p1:-?}, p2_vec_rows=${p2:-?})"; fail=1
	fi
fi
drive "$W/ref" "$H/real_cases.txt" "$W/ref.out"
drive "$W/simd" "$H/real_cases.txt" "$W/simd.out"
same "yama fuzz real cases" "$W/ref.out" "$W/simd.out"
# the run-time override: the dispatch build with YAMA_SIMD=0 must equal the scalar reference too
python3 "$H/gen.py" 7 300 200 > "$W/in"
drive "$W/ref" "$W/in" "$W/ref.out"
drive "$W/simd" "$W/in" "$W/ov.out" YAMA_SIMD=0
same "yama fuzz override YAMA_SIMD=0" "$W/ref.out" "$W/ov.out"
exit $fail
