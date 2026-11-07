#!/bin/sh
# Regression tests: tests/run.sh [-g] [bindir]
# Runs the programs in bindir (default: the source directory) on the inputs
# in tests/data and compares the output, '#' lines aside, with
# tests/expected. -g writes tests/expected instead of comparing.
# Prints one PASS or FAIL line per case and a tally; exits 1 on any FAIL.

T=$(cd "$(dirname "$0")" && pwd)
gen=0
[ "${1-}" = -g ] && { gen=1; shift; }
B=$(cd "${1:-$T/..}" && pwd) || exit 2
D=$T/data
E=$T/expected
W=$(mktemp -d "${TMPDIR:-/tmp}/multiz-check.XXXXXX") || exit 2
trap 'rm -rf "$W"' EXIT
trap 'exit 2' HUP INT TERM
PATH=$B:$PATH
export PATH
TREE="((sa sb) (sc (sd se)))"
OFFSETS="2147480000 5000000000 4000000000000"
pass=0
fail=0

# result NAME WHY: WHY is empty on success
result() {
	if [ -z "$2" ]; then
		echo "PASS $1"
		pass=$((pass + 1))
	else
		echo "FAIL $1:$2"
		fail=$((fail + 1))
	fi
}

# same GOT EXPECTED: equal once '#' lines are dropped; -g copies GOT
same() {
	if [ $gen = 1 ]; then
		grep -v '^#' "$1" > "$2"
		return 0
	fi
	grep -v '^#' "$1" > "$W/got"
	cmp -s "$W/got" "$2"
}

# shift OFF IN OUT
shift_maf() {
	awk -v off="$1" -f "$T/shift.awk" "$2" > "$3"
}

# multiz DIR A B: merge into DIR/out.maf, DIR/U1.maf, DIR/U2.maf; DIR/rc
multiz_run() {
	mkdir -p "$1"
	(cd "$1" && multiz M=1 "$2" "$3" 0 U1.maf U2.maf > out.maf 2> err; echo $? > rc)
}

# roast DIR: DIR holds sa.*.sing.maf
roast_run() {
	mkdir -p "$1/tmp"
	(cd "$1" && roast T="$1/tmp" E=sa "$TREE" sa.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
}

# ran DIR: exited 0 and wrote at least one block
ran() {
	[ "$(cat "$1/rc")" = 0 ] && grep -q '^a' "$1/out.maf"
}

# 1. multiz on synthetic pairs, both orders, against the expected output.
for p in k2l2s1 k8l8s2 k3l40s1; do
	for o in fwd rev; do
		a=$D/synth/$p.1.maf b=$D/synth/$p.2.maf
		[ $o = rev ] && { a=$D/synth/$p.2.maf b=$D/synth/$p.1.maf; }
		d=$W/$p.$o
		multiz_run "$d" "$a" "$b"
		why=
		ran "$d" || why=" exit $(cat "$d/rc")"
		for f in out U1 U2; do
			same "$d/$f.maf" "$E/$p.$o.$f.maf" || why="$why $f.maf differs"
		done
		result "multiz $p $o" "$why"
	done
done

# 2. roast, then tba, on lastz alignments of five synthetic species.
d=$W/roast
mkdir -p "$d"
cp "$D"/tba/sa.*.sing.maf "$d"
roast_run "$d"
why=
ran "$d" || why=" exit $(cat "$d/rc")"
same "$d/out.maf" "$E/roast.maf" || why="$why out.maf differs"
result "roast" "$why"

d=$W/tba
mkdir -p "$d/tmp"
cp "$D"/tba/* "$d"
(cd "$d" && tba T="$d/tmp" "$TREE" *.*.sing.maf tba.maf > /dev/null 2> err; echo $? > rc)
why=
[ "$(cat "$d/rc")" = 0 ] && grep -q '^a' "$d/tba.maf" || why=" exit $(cat "$d/rc")"
same "$d/tba.maf" "$E/tba.maf" || why="$why tba.maf differs"
result "tba" "$why"

[ $gen = 1 ] && { echo "wrote $E"; exit 0; }

# 3. Shift invariance: every coordinate moved past 2^31 and 2^32 must give
# the unshifted output, moved. Needs no expected files.
for off in $OFFSETS; do
	for p in k2l2s1 k8l8s2 k3l40s1; do
		d=$W/$p.s$off
		mkdir -p "$d"
		shift_maf "$off" "$D/synth/$p.1.maf" "$d/1.maf"
		shift_maf "$off" "$D/synth/$p.2.maf" "$d/2.maf"
		multiz_run "$d/run" "$d/1.maf" "$d/2.maf"
		why=
		ran "$d/run" || why=" exit $(cat "$d/run/rc")"
		for f in out U1 U2; do
			shift_maf "$off" "$W/$p.fwd/$f.maf" "$d/want"
			shift_maf 0 "$d/run/$f.maf" "$d/got"
			cmp -s "$d/want" "$d/got" || why="$why $f.maf differs"
		done
		result "multiz $p shifted $off" "$why"
	done
	d=$W/roast.s$off
	mkdir -p "$d"
	for f in "$D"/tba/sa.*.sing.maf; do
		shift_maf "$off" "$f" "$d/${f##*/}"
	done
	roast_run "$d"
	why=
	ran "$d" || why=" exit $(cat "$d/rc")"
	shift_maf "$off" "$W/roast/out.maf" "$d/want"
	shift_maf 0 "$d/out.maf" "$d/got"
	cmp -s "$d/want" "$d/got" || why="$why out.maf differs"
	result "roast shifted $off" "$why"
done

echo "$pass passed, $fail failed"
[ $fail = 0 ]
