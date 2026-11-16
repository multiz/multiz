#!/usr/bin/env bash
# Mutation test of the AVX2 kernels in mz_yama.c.
#   tests/simd_mutation/run_mutants.sh WORKDIR BASE_TREE [-j N] [--no-piece] [MUTANT...]
# BASE_TREE is a scalar multiz tree (e.g. `git archive v11.9`), the fuzzer's reference.
# Mutants are listed in mutants.py (default: all).  The unmutated tree, "none", runs first
# as the control: if any test fails on it, the run is void and the script exits 2.
# For each mutant: copy this tree, apply the mutant (mutants.py exits 1 unless it applied
# exactly once), build with the Makefile's flags, then:
#   runsh   tests/run.sh
#   fuzz    tests/yama_fuzz/run.sh against BASE_TREE: seeds 1-8 x 400 cases; seeds 301-304
#           with YAMA_FUZZ_GOP=268435456 (row coefficients past 32 bits); and seeds 401-404
#           with YAMA_FUZZ_GOP_IDX=2, 0, then 8, which enlarge only Ix, only Cx, or only
#           Cx and Cz (see tests/yama_fuzz/main.c), so no FITS32 check masks another
#   shadow  the same fuzz runs built -DYAMA_SHADOW; killed on a cell MISMATCH, an output
#           DIFF, or zero rows through either vector kernel
#   piece   hg38 100-way piece 03, '#' lines dropped, against its digest (--no-piece skips)
# Every check runs under timeout(1): TEST_TIMEOUT seconds (default 900) for runsh, each
# fuzz run and each shadow run, PIECE_TIMEOUT (default 5400) for piece.  A timeout is a
# kill (a mutant that hangs is detected) and is named in the notes.
# Writes WORKDIR/results.tsv: mutant runsh fuzz shadow piece notes, with K (killed),
# S (survived), B (did not build) or - (not run).  Exits 1 if any mutant survived every test.
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$H/../.." && pwd)
FZ=$SRC/tests/yama_fuzz/run.sh
PIECE=${PIECE:-/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts/ucsc100way_run_piece.sh}
U=${U:-/home/rcb112/project/claude/bench/data/ucsc100way}
DIGEST=de0dca97a159c7a0587e29f67a3537a77427357eb523f7882c03c5460d4374db
GOP=268435456
TO=${TEST_TIMEOUT:-900}
PTO=${PIECE_TIMEOUT:-5400}

if [ "${1-}" = --one ]; then	# internal: one mutant
	m=$2 W=$3 BT=$4 dopiece=$5
	d=$W/m_$m
	rm -rf "$d"; mkdir -p "$d/src" "$d/tmp"
	(cd "$SRC" && git ls-files | tar -cf - -T -) | tar -xf - -C "$d/src"
	notes=""
	if [ "$m" != none ]; then
		python3 -I "$H/mutants.py" apply "$m" "$d/src/mz_yama.c" || { printf '%s\t?\t?\t?\t?\tmutant did not apply\n' "$m" > "$d/result"; exit 0; }
	fi
	if ! (cd "$d/src" && make -j2 > "$d/build.log" 2>&1); then
		printf '%s\tB\tB\tB\tB\tbuild failed\n' "$m" > "$d/result"; exit 0
	fi
	# chk TEST SECONDS CMD...: run CMD under timeout; on failure mark TEST killed (eval'd
	# into the variable of that name), and note a timeout (exit 124 or 137)
	chk() {
		local t=$1 secs=$2 rc; shift 2
		timeout -k 10 "$secs" "$@"; rc=$?
		[ $rc = 0 ] && return 0
		eval "$t=K"
		case $rc in 124|137) notes="$notes timeout:$t" ;; esac
		return 1
	}
	runsh=S fuzz=S shadow=S
	chk runsh "$TO" env TMPDIR="$d/tmp" sh "$SRC/tests/run.sh" "$d/src" > "$d/runsh.out" 2>&1
	chk fuzz "$TO" sh "$FZ" "$BT" "$d/src" "$d/fz" "1 2 3 4 5 6 7 8" 400 200 -O2 > "$d/fuzz.out" 2>&1
	chk fuzz "$TO" env YAMA_FUZZ_GOP=$GOP sh "$FZ" "$BT" "$d/src" "$d/fzg" "301 302 303 304" 400 200 -O2 >> "$d/fuzz.out" 2>&1
	for ix in 2 0 8; do
		chk fuzz "$TO" env YAMA_FUZZ_GOP=$GOP YAMA_FUZZ_GOP_IDX=$ix sh "$FZ" "$BT" "$d/src" "$d/fzi$ix" "401 402 403 404" 400 200 -O2 >> "$d/fuzz.out" 2>&1
	done
	rm -f "$d/shadow.log"
	export YAMA_SHADOW_LOG=$d/shadow.log
	chk shadow "$TO" sh "$FZ" "$BT" "$d/src" "$d/sfz" "1 2 3 4 5 6 7 8" 400 200 "-O2 -DYAMA_SHADOW" > "$d/shadow.out" 2>&1
	chk shadow "$TO" env YAMA_FUZZ_GOP=$GOP sh "$FZ" "$BT" "$d/src" "$d/sfzg" "301 302 303 304" 400 200 "-O2 -DYAMA_SHADOW" >> "$d/shadow.out" 2>&1
	for ix in 2 0 8; do
		chk shadow "$TO" env YAMA_FUZZ_GOP=$GOP YAMA_FUZZ_GOP_IDX=$ix sh "$FZ" "$BT" "$d/src" "$d/sfzi$ix" "401 402 403 404" 400 200 "-O2 -DYAMA_SHADOW" >> "$d/shadow.out" 2>&1
	done
	unset YAMA_SHADOW_LOG
	if grep -q MISMATCH "$d/shadow.log" 2>/dev/null; then
		shadow=K; notes="$notes shadow:$(grep -m1 -o 'field=[A-Za-z_]*' "$d/shadow.log")"
	fi
	vr=$(awk '{for(i=2;i<=NF;i++){split($i,a,"=");if(a[1]=="p1_vec_rows")p+=a[2];if(a[1]=="p2_vec_rows")q+=a[2]}} END{print (p>0&&q>0)?"ok":"zero"}' "$d/shadow.log" 2>/dev/null)
	[ "$vr" = ok ] || { shadow=K; notes="$notes shadow:no-vector-rows"; }
	piece=-
	if [ "$dopiece" = 1 ]; then
		piece=S
		if chk piece "$PTO" bash "$PIECE" "$d/src" "$U/run03/inputs-v3" "$U/tree.nh" "$d/piece" > "$d/piece.out" 2>&1; then
			[ "$(grep -v '^#' "$d/piece/out.maf" | sha256sum | cut -c1-64)" = $DIGEST ] || piece=K
		fi
		rm -rf "$d/piece/tmp" "$d/piece/out.maf"
	fi
	rm -rf "$d"/fz* "$d"/sfz* "$d/tmp"
	printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$runsh" "$fuzz" "$shadow" "$piece" "${notes# }" > "$d/result"
	exit 0
fi

W=$1 BT=$(cd "$2" && pwd); shift 2
J=4 dopiece=1 names=()
while [ $# -gt 0 ]; do
	case $1 in
	-j) J=$2; shift 2 ;;
	--no-piece) dopiece=0; shift ;;
	*) names+=("$1"); shift ;;
	esac
done
[ ${#names[@]} -gt 0 ] || mapfile -t names < <(python3 -I "$H/mutants.py" list)
mkdir -p "$W"; W=$(cd "$W" && pwd)

"$0" --one none "$W" "$BT" "$dopiece"
IFS=$'\t' read -r _ a b c e _ < "$W/m_none/result"
if [ "$a$b$c" != SSS ] || { [ "$dopiece" = 1 ] && [ "$e" != S ]; }; then
	echo "control (unmutated tree) failed: $(cat "$W/m_none/result")"; exit 2
fi
printf '%s\n' "${names[@]}" | xargs -P "$J" -I{} "$0" --one {} "$W" "$BT" "$dopiece"
{
	printf 'mutant\trunsh\tfuzz\tshadow\tpiece\tnotes\n'
	cat "$W/m_none/result"
	for m in "${names[@]}"; do cat "$W/m_$m/result"; done
} > "$W/results.tsv"
column -t -s $'\t' "$W/results.tsv"
awk -F'\t' 'NR>2 && $2!="K" && $3!="K" && $4!="K" && $5!="K" && $2!="B" {print "SURVIVOR", $1; s=1} END{exit s}' "$W/results.tsv"
