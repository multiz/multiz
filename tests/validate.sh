#!/usr/bin/env bash
# Release validation of the AVX2 yama: sanitizers, valgrind, and a compiler/optimization matrix.
#
#   tests/validate.sh build NAME CC "CFLAGS"     build this tree (HEAD plus the working tree) into $VW/NAME
#   tests/validate.sh suite NAME [ENV...]        tests/run.sh and tests/yama_fuzz/check.sh against NAME, under ENV
#   tests/validate.sh deep NAME [ENV...]         the cannabis deep merge; prints the digest of its three outputs
#   tests/validate.sh fuzzdig NAME [ENV...]      digest of the yama fuzz driver's output on fixed seeds (built with NAME's CC/CFLAGS)
#   tests/validate.sh valgrind NAME              memcheck on tests/run.sh (every program wrapped) and on the deep merge
#   tests/validate.sh matrix                     every CC x CFLAGS build: suite, fuzz digest, deep digest, all equal
#   tests/validate.sh sanitize                   gcc and clang ASan+UBSan builds: suite and deep merge, both dispatch paths
#   tests/validate.sh sanitize-recover           the same, logging every report and continuing; zero reports from mz_yama.c required
#   tests/validate.sh bounds                     the z3 bounds asserted (-DYAMA_ASSERT_BOUNDS) over every workload; PIECE03=1 adds piece 03
#   tests/validate.sh cross                      static aarch64 build under qemu-aarch64: scalar fallback, no warnings, x86-identical output
#   tests/validate.sh conda                      bioconda's typical CFLAGS (-march=nocona ...): identical output, AVX2 engages at run time
#
# VW (default bench/data/tmp/simd_validate) holds the builds and logs. A sanitizer report makes the program exit
# nonzero (halt_on_error, -fno-sanitize-recover), so a clean suite means zero reports. Leak detection is OFF under
# ASan: multiz has always freed little at exit; `valgrind` compares leak totals of the scalar and AVX2 paths instead.
# The piece-03 workload is too big for this host under ASan; bench/laila_release.sh covers it on laila.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
T=$(cd "$H/.." && pwd)
VW=${VW:-/home/rcb112/project/claude/bench/data/tmp/simd_validate}
DEEP=${DEEP:-/home/rcb112/project/claude/bench/data/wfi_profile/in}
PROGS="all_bz lav2maf pair2tb single_cov2 blastzWrapper maf_project maf2lav maf2fasta maf_order mafFind get_standard_headers maf_checkThread tba roast multic multiz maf_sort get_covered"
export TMPDIR=$VW/tmp
mkdir -p "$TMPDIR"
export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1:abort_on_error=0:exitcode=86}
export UBSAN_OPTIONS=${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1:exitcode=87}
SAN="-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer -g"
# the matrix: name, compiler, flags
MATRIX=(
	"gcc-O0 gcc -O0" "gcc-O2 gcc -O2" "gcc-O3 gcc -O3" "gcc-O3n gcc -O3 -march=native"
	"clang-O0 clang -O0" "clang-O2 clang -O2" "clang-O3 clang -O3" "clang-O3n clang -O3 -march=native"
)

cfg() { cat "$VW/$1/.cfg"; } # "CC|CFLAGS"

build() {
	local name=$1 cc=$2 cf=$3 d=$VW/$1
	rm -rf "$d"; mkdir -p "$d"
	# the tree as it is now: tracked files at their working-tree state
	(cd "$T" && git ls-files -z | xargs -0 tar -cf - --no-recursion) | tar -C "$d" -xf -
	make -s -C "$d" -j8 CC="$cc" CFLAGS="-Wall -Wextra $cf" $PROGS > "$d/.build.log" 2>&1 \
		|| { echo "BUILD FAILED $name"; tail -20 "$d/.build.log"; exit 1; }
	printf '%s|%s\n' "$cc" "$cf" > "$d/.cfg"
	echo "BUILT $name ($cc $cf)"
}

suite() {
	local name=$1; shift
	local c; c=$(cfg "$name")
	local r1=0 r2=0
	env "$@" sh "$T/tests/run.sh" "$VW/$name" > "$VW/$name/.suite.log" 2>&1 || r1=$?
	env "$@" YAMA_FUZZ_CC="${c%%|*}" YAMA_FUZZ_CFLAGS="${c#*|}" sh "$T/tests/yama_fuzz/check.sh" "$VW/$name" \
		> "$VW/$name/.fuzz.log" 2>&1 || r2=$?
	echo "SUITE $name $* run.sh=$(tail -1 "$VW/$name/.suite.log") rc=$r1 fuzz_rc=$r2 fuzz_fail=$(grep -c '^FAIL' "$VW/$name/.fuzz.log" || true)"
	[ $r1 = 0 ] && { [ $r2 = 0 ] || [ $r2 = 77 ]; }
}

deep() {
	local name=$1; shift
	local o=$VW/$name/.deep; rm -rf "$o"; mkdir -p "$o"
	env "$@" "$VW/$name/multiz" M=1 "$DEEP/left.maf" "$DEEP/right.maf" 0 "$o/U1.maf" "$o/U2.maf" > "$o/out.maf" 2> "$o/err" \
		|| { echo "DEEP FAILED $name rc=$?"; tail -20 "$o/err"; return 1; }
	echo "DEEP $name $* $(cat "$o/out.maf" "$o/U1.maf" "$o/U2.maf" | grep -v '^#' | sha256sum | cut -c1-16)"
	rm -f "$o"/*.maf
}

fuzzdig() {
	local name=$1; shift
	local c; c=$(cfg "$name")
	local cc=${c%%|*} cf=${c#*|} f=$VW/$name/.yama_fuzz
	$cc $cf -Wno-free-nonheap-object -I"$VW/$name" "$T/tests/yama_fuzz/main.c" "$VW/$name/mz_yama.c" \
		"$VW/$name/mz_scores.c" "$VW/$name/util.c" -o "$f" 2> /dev/null
	local dig
	dig=$(for s in 11 12 13; do python3 "$T/tests/yama_fuzz/gen.py" $s 300 200 | env "$@" "$f"; done | sha256sum | cut -c1-16)
	echo "FUZZ $name $* $dig"
}

valgrind_run() {
	local name=$1 v=$VW/$1.vg
	rm -rf "$v"; mkdir -p "$v/bin" "$v/logs"
	for p in $PROGS; do
		printf '#!/bin/sh\nexec valgrind -q --error-exitcode=99 --leak-check=no --log-file=%s/logs/%s.%%p "%s/%s" "$@"\n' \
			"$v" "$p" "$VW/$name" "$p" > "$v/bin/$p"
		chmod +x "$v/bin/$p"
	done
	sh "$T/tests/run.sh" "$v/bin" > "$v/suite.log" 2>&1 || true
	local errs; errs=$(cat "$v"/logs/* 2> /dev/null | grep -c 'ERROR SUMMARY\|Invalid\|uninitialised' || true)
	echo "VALGRIND suite $name: $(tail -1 "$v/suite.log"); memcheck error lines: $errs"
	for mode in simd scalar; do
		local envs=(); [ $mode = scalar ] && envs=(YAMA_SIMD=0)
		local rc=0
		env "${envs[@]}" valgrind --error-exitcode=99 --leak-check=summary --log-file="$v/deep.$mode.log" \
			"$VW/$name/multiz" M=1 "$DEEP/left.maf" "$DEEP/right.maf" 0 "$v/U1.maf" "$v/U2.maf" > "$v/out.maf" 2> /dev/null || rc=$?
		echo "VALGRIND deep $name $mode rc=$rc $(grep 'ERROR SUMMARY' "$v/deep.$mode.log" | sed 's/==[0-9]*==//') $(grep -E 'definitely lost|indirectly lost' "$v/deep.$mode.log" | sed 's/==[0-9]*==//' | tr -s ' ' | paste -sd';')"
		rm -f "$v"/*.maf
	done
}

cmd=${1:?usage: see the header}; shift
case $cmd in
build) build "$@" ;;
suite) suite "$@" ;;
deep) deep "$@" ;;
fuzzdig) fuzzdig "$@" ;;
valgrind) valgrind_run "$@" ;;
matrix)
	rc=0
	for m in "${MATRIX[@]}"; do set -- $m; n=$1 cc=$2; shift 2; build "$n" "$cc" "$*"; done
	for m in "${MATRIX[@]}"; do
		n=${m%% *}
		suite "$n" || rc=1
		fuzzdig "$n"; fuzzdig "$n" YAMA_SIMD=0
		deep "$n"
	done | tee "$VW/matrix.log"
	# every build and both paths must agree on the fuzz digest, and every build on the deep digest
	for k in FUZZ DEEP; do
		u=$(awk -v k=$k '$1 == k {print $NF}' "$VW/matrix.log" | sort -u | wc -l)
		echo "MATRIX $k distinct digests: $u"; [ "$u" = 1 ] || rc=1
	done
	grep -q '^SUITE.*rc=0 fuzz_rc=0' "$VW/matrix.log" || rc=1
	[ "$(grep -c '^SUITE' "$VW/matrix.log")" = "$(grep -c '^SUITE.*rc=0 fuzz_rc=0' "$VW/matrix.log")" ] || rc=1
	exit $rc
	;;
sanitize)
	rc=0
	for c in gcc clang; do
		build "san-$c" "$c" "-O1 $SAN"
		suite "san-$c" || rc=1
		suite "san-$c" YAMA_SIMD=0 || rc=1
		deep "san-$c" || rc=1
		deep "san-$c" YAMA_SIMD=0 || rc=1
	done 2>&1 | tee "$VW/sanitize.log"
	grep -q 'FAILED\|rc=[1-9]\|fuzz_rc=[1-9]' "$VW/sanitize.log" && rc=1
	exit $rc
	;;
sanitize-recover)
	# Like `sanitize`, but every report is logged and execution continues, so one pre-existing error elsewhere
	# (mz_preyama.c reads map4[-1]; see simd_row_prototype.md) cannot hide later ones. Each report is attributed to
	# its first frame in this tree; the requirement is ZERO reports from mz_yama.c, on both dispatch paths.
	rc=0
	command -v llvm-symbolizer > /dev/null && export ASAN_SYMBOLIZER_PATH=$(command -v llvm-symbolizer)
	for c in gcc clang; do
		n=sanr-$c
		build "$n" "$c" "-O1 -fsanitize=address,undefined -fsanitize-recover=all -fno-omit-frame-pointer -g"
		for mode in simd scalar; do
			envs=(); [ $mode = scalar ] && envs=(YAMA_SIMD=0)
			L=$VW/$n/.san.$mode; rm -rf "$L"; mkdir -p "$L"
			export ASAN_OPTIONS="detect_leaks=0:halt_on_error=0:log_path=$L/asan"
			export UBSAN_OPTIONS="halt_on_error=0:print_stacktrace=1:log_path=$L/ubsan"
			suite "$n" "${envs[@]}" || true
			deep "$n" "${envs[@]}" || true
			# tally: the first "#k ... in FUNC FILE:LINE" frame whose FILE is in the build tree, per report
			# with no report at all there are no log files, and awk must not then read stdin
			reps=(); for f in "$L"/asan* "$L"/ubsan*; do [ -e "$f" ] && reps+=("$f"); done
			[ ${#reps[@]} -gt 0 ] || reps=(/dev/null)
			awk -v root="$VW/$n/" '
				/ERROR: AddressSanitizer|runtime error/ { want = 1; kind = ($0 ~ /runtime error/) ? "ubsan" : "asan"
					if (kind == "ubsan") { split($1, a, ":"); f = a[1]; sub(root, "", f); if (f != "") { print kind, f ":" a[2]; want = 0 } } next }
				want && /^ *#[0-9]+ / && index($0, root) { for (i = 1; i <= NF; i++) if (index($i, root)) { f = $i; sub(root, "", f); sub(/:[0-9]+$/, "", f); print kind, f; want = 0; break } }
			' "${reps[@]}" | sort | uniq -c | sort -rn > "$L/tally"
			echo "TALLY $n $mode: $(awk '{printf "%s %s x%d; ", $2, $3, $1}' "$L/tally")"
			if grep -q 'mz_yama\.c' "$L/tally"; then echo "FAIL $n $mode: a report from mz_yama.c"; rc=1; fi
		done
	done 2>&1 | tee "$VW/sanitize_recover.log"
	grep -q '^FAIL' "$VW/sanitize_recover.log" && rc=1
	exit $rc
	;;
bounds)
	# The z3 proof's value bounds, asserted (-DYAMA_ASSERT_BOUNDS, which also runs shadow mode) over every workload:
	# tests/run.sh, the deep merge, every check.sh input set (random, 2^28 gap opens, single-entry tables, targeted,
	# real calls) and, with PIECE03=1, hg38 100-way piece 03. Alignments (stdout) must equal a plain scalar build;
	# shadow mode's counters go to stderr and differ by design. Prints the observed extremes.
	rc=0; F=$T/tests/yama_fuzz; L=$VW/bounds.shadow.log; rm -f "$L"
	build ab gcc "-O2 -DYAMA_ASSERT_BOUNDS"
	export YAMA_SHADOW_LOG=$L
	gcc -O2 -Wno-free-nonheap-object -DYAMA_ASSERT_BOUNDS -I"$T" "$F/main.c" "$T/mz_yama.c" "$T/mz_scores.c" "$T/util.c" -o "$VW/ab.fz"
	gcc -O2 -Wno-free-nonheap-object -DYAMA_SIMD=0 -I"$T" "$F/main.c" "$T/mz_yama.c" "$T/mz_scores.c" "$T/util.c" -o "$VW/ref.fz"
	r=0; sh "$T/tests/run.sh" "$VW/ab" > "$VW/ab/.suite.log" 2>&1 || r=$?
	echo "BOUNDS run.sh: $(tail -1 "$VW/ab/.suite.log") rc=$r"; [ $r = 0 ] || rc=1
	deep ab || rc=1
	fz() { # LABEL GENARGS... -- ENV...
		local label=$1; shift; local ga=(); while [ "$1" != -- ]; do ga+=("$1"); shift; done; shift
		if [ "${ga[0]}" = real ]; then cp "$F/real_cases.txt" "$VW/b.in"; else python3 "$F/gen.py" "${ga[@]}" > "$VW/b.in"; fi
		local a=0; env "$@" "$VW/ab.fz" < "$VW/b.in" > "$VW/b.ab" 2> /dev/null || a=$?
		env "$@" "$VW/ref.fz" < "$VW/b.in" > "$VW/b.ref" 2> /dev/null
		if [ $a = 0 ] && cmp -s "$VW/b.ab" "$VW/b.ref"; then echo "BOUNDS fuzz $label: ok"; else echo "BOUNDS fuzz $label: FAIL rc=$a"; rc=1; fi
	}
	G=268435456
	for s in 1 2 3; do fz "random $s" $s 300 200 --; done
	for s in 301 302; do fz "gop2^28 $s" $s 300 200 -- YAMA_FUZZ_GOP=$G; done
	for ix in 2 0 8; do fz "gop_idx $ix" 401 300 200 -- YAMA_FUZZ_GOP=$G YAMA_FUZZ_GOP_IDX=$ix; done
	for m in narrow ties edges; do fz "targeted $m" 501 300 200 $m --; done
	fz "real cases" real --
	if [ "${PIECE03:-0}" = 1 ]; then
		P=${PIECE:-/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts/ucsc100way_run_piece.sh}
		U=${U:-/home/rcb112/project/claude/bench/data/ucsc100way}
		rm -rf "$VW/ab.piece"; bash "$P" "$VW/ab" "$U/run03/inputs-v3" "$U/tree.nh" "$VW/ab.piece" > /dev/null 2>&1 || true
		dg=$(grep -v '^#' "$VW/ab.piece/out.maf" 2> /dev/null | sha256sum | cut -c1-64)
		echo "BOUNDS piece03 digest $dg (expected de0dca97a159c7a0587e29f67a3537a77427357eb523f7882c03c5460d4374db)"
		[ "$dg" = de0dca97a159c7a0587e29f67a3537a77427357eb523f7882c03c5460d4374db ] || rc=1
	fi
	n=$(grep -c 'YAMA_ASSERT_BOUNDS\|MISMATCH' "$L" 2> /dev/null || true)
	awk '/^YAMA_SHADOW calls/ { for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
		if (!seen || v["vmin"] + 0 < vmin) vmin = v["vmin"] + 0; if (!seen || v["vmax"] + 0 > vmax) vmax = v["vmax"] + 0
		if (v["vmin_above_MININT"] + 0 < vmn || !seen) vmn = v["vmin_above_MININT"] + 0
		if (v["max_abs_Izc"] + 0 > z) z = v["max_abs_Izc"] + 0; if (v["max_abs_eK"] + 0 > e) e = v["max_abs_eK"] + 0
		cells += v["cells"]; seen = 1 }
		END { printf "BOUNDS extremes over %d processes: vmin=%.0f vmax=%.0f vmin_above_MININT=%.0f max|Izc|=%.0f max|eK|=%.0f cells=%.0f\n", NR, vmin, vmax, vmn, z, e, cells }' "$L"
	echo "BOUNDS assertion or mismatch lines: $n"; [ "$n" = 0 ] || rc=1
	exit $rc
	;;
cross)
	# Non-x86 (bioconda builds linux-aarch64 and osx-arm64): a static aarch64 build run under qemu-aarch64.
	# YAMA_SIMD must evaluate to 0 (the AVX2 code compiled out), mz_yama.c must compile warning-free, and the output
	# must equal the x86 builds': tests/run.sh, check.sh, the fuzz digest and the deep-merge digest.
	rc=0; X=aarch64-linux-gnu-gcc; Q=qemu-aarch64; n=x-aarch64
	build $n $X "-O2 -static"
	w=$(grep -c 'mz_yama\.c.*warning' "$VW/$n/.build.log" || true)
	v=$(echo '#include "mz_yama.c"' | $X -O2 -dM -E -I"$T" -x c - 2> /dev/null | awk '$2 == "YAMA_SIMD" {print $3}')
	echo "CROSS $n: mz_yama.c warnings $w; YAMA_SIMD=$v"; [ "$w" = 0 ] && [ "$v" = 0 ] || rc=1
	b=$VW/$n.q; rm -rf "$b"; mkdir -p "$b"
	for p in $PROGS; do printf '#!/bin/sh\nexec %s "%s/%s" "$@"\n' "$Q" "$VW/$n" "$p" > "$b/$p"; chmod +x "$b/$p"; done
	r=0; sh "$T/tests/run.sh" "$b" > "$VW/$n/.suite.log" 2>&1 || r=$?
	echo "CROSS run.sh: $(tail -1 "$VW/$n/.suite.log") rc=$r"; [ $r = 0 ] || rc=1
	r=0; YAMA_FUZZ_CC=$X YAMA_FUZZ_CFLAGS="-O2 -static" YAMA_FUZZ_RUN=$Q sh "$T/tests/yama_fuzz/check.sh" "$VW/$n" \
		> "$VW/$n/.fuzz.log" 2>&1 || r=$?
	echo "CROSS check.sh: rc=$r $(grep -c '^PASS' "$VW/$n/.fuzz.log") PASS $(grep -c '^FAIL' "$VW/$n/.fuzz.log") FAIL"; [ $r = 0 ] || rc=1
	f=$VW/$n/.yama_fuzz
	$X -O2 -static -Wno-free-nonheap-object -I"$VW/$n" "$T/tests/yama_fuzz/main.c" "$VW/$n/mz_yama.c" "$VW/$n/mz_scores.c" \
		"$VW/$n/util.c" -o "$f"
	dg=$(for s in 11 12 13; do python3 "$T/tests/yama_fuzz/gen.py" $s 300 200 | $Q "$f"; done | sha256sum | cut -c1-16)
	echo "CROSS fuzz digest $dg (x86 matrix: ${X86_FUZZ:-b3d5542704b57607})"; [ "$dg" = "${X86_FUZZ:-b3d5542704b57607}" ] || rc=1
	o=$VW/$n/.deep; rm -rf "$o"; mkdir -p "$o"
	$Q "$VW/$n/multiz" M=1 "$DEEP/left.maf" "$DEEP/right.maf" 0 "$o/U1.maf" "$o/U2.maf" > "$o/out.maf" 2> "$o/err" || rc=1
	dd=$(cat "$o/out.maf" "$o/U1.maf" "$o/U2.maf" | grep -v '^#' | sha256sum | cut -c1-16)
	echo "CROSS deep digest $dd (x86 matrix: ${X86_DEEP:-066d7c813aa8ed1c})"; [ "$dd" = "${X86_DEEP:-066d7c813aa8ed1c}" ] || rc=1
	rm -f "$o"/*.maf
	exit $rc
	;;
conda)
	# Bioconda's typical x86-64 CFLAGS. -march=nocona has no AVX2, so this checks that the target("avx2") kernels
	# still compile and, on AVX2 hardware, ENGAGE at run time (shadow mode counts the rows each kernel did), and that
	# the output equals the matrix builds'.
	rc=0; CF="-march=nocona -mtune=haswell -O2 -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe"
	build conda gcc "$CF"
	w=$(grep -c 'mz_yama\.c.*warning' "$VW/conda/.build.log" || true); echo "CONDA mz_yama.c warnings: $w"; [ "$w" = 0 ] || rc=1
	suite conda || rc=1
	fuzzdig conda | tee "$VW/conda.fz"; fuzzdig conda YAMA_SIMD=0
	deep conda | tee "$VW/conda.dp"
	grep -q "${X86_FUZZ:-b3d5542704b57607}" "$VW/conda.fz" || rc=1; grep -q "${X86_DEEP:-066d7c813aa8ed1c}" "$VW/conda.dp" || rc=1
	build conda-shadow gcc "$CF -DYAMA_SHADOW"
	L=$VW/conda.shadow.log; rm -f "$L"
	YAMA_SHADOW_LOG=$L "$VW/conda-shadow/multiz" M=1 "$DEEP/left.maf" "$DEEP/right.maf" 0 /dev/null /dev/null > /dev/null 2>&1 || rc=1
	awk '/^YAMA_SHADOW calls/ { for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] += kv[2] } }
		END { printf "CONDA AVX2 engaged on the deep merge: p1_vec_rows=%d p2_vec_rows=%d of rows=%d (cells %d, vector p1 %d)\n", v["p1_vec_rows"], v["p2_vec_rows"], v["rows"], v["cells"], v["p1_vec_cells"] }' "$L" | tee "$VW/conda.eng"
	# p1_vec_rows by name: field 6 is the word "deep", so the old test always failed, hiding the digest check above
	awk '{ for (i = 1; i <= NF; i++) if ($i ~ /^p1_vec_rows=/) { split($i, a, "="); exit !(a[2] > 0) } exit 1 }' "$VW/conda.eng" || rc=1
	grep -q MISMATCH "$L" && rc=1
	exit $rc
	;;
*) echo "unknown command $cmd" >&2; exit 2 ;;
esac
