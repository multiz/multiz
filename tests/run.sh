#!/bin/sh
# Regression tests: tests/run.sh [-g] [bindir]
# Runs the programs in bindir (default: the source directory) on the inputs
# in tests/data and compares the output, '#' lines aside, with
# tests/expected. -g writes tests/expected instead of comparing, and only
# when every case it generates from ran cleanly; otherwise it changes nothing.
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
# the experimental option (section 2b) must not leak into the default cases from the caller's environment
unset MULTIZ_PREYAMA_SKIP_DASH
TREE="((sa sb) (sc (sd se)))"
UTREE="((hg38 allMis1) (galGal4 taeGut2))"
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

# same GOT EXPECTED: equal once '#' lines are dropped; -g stages GOT as EXPECTED.gen, which is moved into
# place only if the whole -g run passes (a failed run must not become the expected output)
same() {
	if [ $gen = 1 ]; then
		grep -v '^#' "$1" > "$2.gen"
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

d=$W/ucsc
mkdir -p "$d/tmp"
cp "$D"/ucsc/*.sing.maf "$d"
(cd "$d" && roast T="$d/tmp" E=hg38 "$UTREE" hg38.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
why=
ran "$d" || why=" exit $(cat "$d/rc")"
same "$d/out.maf" "$E/roast-ucsc.maf" || why="$why out.maf differs"
result "roast ucsc" "$why"

d=$W/tba
mkdir -p "$d/tmp"
cp "$D"/tba/* "$d"
(cd "$d" && tba T="$d/tmp" "$TREE" *.*.sing.maf tba.maf > /dev/null 2> err; echo $? > rc)
why=
[ "$(cat "$d/rc")" = 0 ] && grep -q '^a' "$d/tba.maf" || why=" exit $(cat "$d/rc")"
same "$d/tba.maf" "$E/tba.maf" || why="$why tba.maf differs"
result "tba" "$why"

# 2b. MULTIZ_PREYAMA_SKIP_DASH, the experimental opt-in pre_yama option: roast E=hg38 on data/skipdash, an
# excerpt of UCSC's hg38 alignments with taeGut2 and zonAlb1 in which the option changes one block. Unset and
# empty give the default output, 1 the option's; the two expected files must differ. Any other value (0 here)
# is off: the default output, with a warning on stderr (before 11.10's review it was fatal, and roast exited 0
# with an empty alignment).
STREE="(hg38 (zonAlb1 taeGut2))"
for m in default empty zero on; do
	d=$W/skipdash.$m
	mkdir -p "$d/tmp"
	cp "$D"/skipdash/*.sing.maf "$d"
	(
		cd "$d" || exit
		case $m in empty) MULTIZ_PREYAMA_SKIP_DASH=; export MULTIZ_PREYAMA_SKIP_DASH ;;
		zero) MULTIZ_PREYAMA_SKIP_DASH=0; export MULTIZ_PREYAMA_SKIP_DASH ;;
		on) MULTIZ_PREYAMA_SKIP_DASH=1; export MULTIZ_PREYAMA_SKIP_DASH ;; esac
		roast T="$d/tmp" E=hg38 "$STREE" hg38.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc
	)
	exp=$E/roast-skipdash.default.maf
	[ $m = on ] && exp=$E/roast-skipdash.maf
	why=
	ran "$d" || why=" exit $(cat "$d/rc")"
	same "$d/out.maf" "$exp" || why="$why out.maf differs"
	if [ $m = on ]; then
		grep -q 'MULTIZ_PREYAMA_SKIP_DASH=1' "$d/err" || why="$why no stderr notice"
		[ $gen = 1 ] || ! cmp -s "$E/roast-skipdash.default.maf" "$E/roast-skipdash.maf" || why="$why the option changed nothing"
	elif [ $m = zero ]; then
		grep -q 'MULTIZ_PREYAMA_SKIP_DASH="0" is not 1' "$d/err" || why="$why no warning"
		! grep -q 'MULTIZ_PREYAMA_SKIP_DASH=1' "$d/err" || why="$why option notice for 0"
	else
		! grep -q 'MULTIZ_PREYAMA_SKIP_DASH' "$d/err" || why="$why notice without the option"
	fi
	result "roast skipdash $m" "$why"
done

if [ $gen = 1 ]; then
	if [ $fail != 0 ]; then
		rm -f "$E"/*.gen
		echo "$fail case(s) failed; $E left unchanged"
		exit 1
	fi
	for f in "$E"/*.gen; do mv "$f" "${f%.gen}"; done
	echo "wrote $E"
	exit 0
fi

# 2c. MULTIZ_PREYAMA_SKIP_DASH in multiz itself: with 1, a '#' line in the header says the option is on (the
# default output has none); any value other than unset, empty or 1 is off: exit 0, the default output, and a
# warning on stderr.
d=$W/skipdash.multiz
(MULTIZ_PREYAMA_SKIP_DASH=1; export MULTIZ_PREYAMA_SKIP_DASH; multiz_run "$d" "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf")
why=
ran "$d" || why=" exit $(cat "$d/rc")"
grep -q '^# MULTIZ_PREYAMA_SKIP_DASH=1' "$d/out.maf" || why="$why no header line"
! grep -q 'MULTIZ_PREYAMA_SKIP_DASH' "$W/k2l2s1.fwd/out.maf" || why="$why header line in the default output"
result "multiz skipdash header" "$why"
for bad in yes 2 0; do
	d=$W/skipdash.bad$bad
	(MULTIZ_PREYAMA_SKIP_DASH=$bad; export MULTIZ_PREYAMA_SKIP_DASH; multiz_run "$d" "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf")
	why=
	ran "$d" || why=" exit $(cat "$d/rc")"
	grep -q "MULTIZ_PREYAMA_SKIP_DASH=\"$bad\" is not 1" "$d/err" || why="$why no warning"
	for f in out U1 U2; do
		grep -v '^#' "$d/$f.maf" > "$W/got"
		grep -v '^#' "$W/k2l2s1.fwd/$f.maf" | cmp -s "$W/got" - || why="$why $f.maf differs from the default"
	done
	result "multiz skipdash $bad is off" "$why"
done

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

# 4. Sequence names containing "eof" or "maf". roast filters those strings out
# of comment lines; an alignment row whose name merely contains one must be
# kept. seofb exercises the merge of the reference with one leaf, where roast
# appends that leaf's pairwise file through the eof filter. Renamed back, the
# output must equal the roast output above. Needs no expected files.
for nm in seofb smafb; do
	d=$W/roast.$nm
	mkdir -p "$d/tmp"
	for f in "$D"/tba/sa.*.sing.maf; do
		g=$(echo "${f##*/}" | sed "s/\.sb\./.$nm./")
		sed "s/^s sb\([ .]\)/s $nm\1/" "$f" > "$d/$g"
	done
	t=$(echo "$TREE" | sed "s/sb/$nm/")
	(cd "$d" && roast T="$d/tmp" E=sa "$t" sa.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
	why=
	ran "$d" || why=" exit $(cat "$d/rc")"
	grep -q "^s $nm" "$d/out.maf" || why="$why no $nm rows in the output"
	# names pad the row columns, so compare fields, not spacing
	sed "s/^s $nm\([ .]\)/s sb\1/" "$d/out.maf" | grep -v '^#' | awk '{$1 = $1; print}' > "$d/got"
	grep -v '^#' "$W/roast/out.maf" | awk '{$1 = $1; print}' > "$d/want"
	cmp -s "$d/want" "$d/got" || why="$why rows lost or changed"
	result "roast name $nm" "$why"
done

# The same on the UCSC excerpt, whose merges leave unused blocks, so the maf
# filter sees alignment rows: allMis1 renamed amafMis1.
nm=amafMis1
d=$W/ucsc.$nm
mkdir -p "$d/tmp"
for f in "$D"/ucsc/*.sing.maf; do
	g=$(echo "${f##*/}" | sed "s/\.allMis1\./.$nm./")
	sed "s/^s allMis1\([ .]\)/s $nm\1/" "$f" > "$d/$g"
done
t=$(echo "$UTREE" | sed "s/allMis1/$nm/")
(cd "$d" && roast T="$d/tmp" E=hg38 "$t" hg38.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
why=
ran "$d" || why=" exit $(cat "$d/rc")"
grep -q "^s $nm" "$d/out.maf" || why="$why no $nm rows in the output"
sed "s/^s $nm\([ .]\)/s allMis1\1/" "$d/out.maf" | grep -v '^#' | awk '{$1 = $1; print}' > "$d/got"
grep -v '^#' "$W/ucsc/out.maf" | awk '{$1 = $1; print}' > "$d/want"
cmp -s "$d/want" "$d/got" || why="$why rows lost or changed"
result "roast ucsc name $nm" "$why"

# 5. Failures are failures (11.10 code review): each case below exited 0, crashed or overran a buffer before.
# None changes the output on valid input.

# a failed write or an unopenable out1 is an error, not exit 0 (finding 4)
d=$W/fail.full
mkdir -p "$d"
(multiz "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf" 0 > /dev/full 2> "$d/err"; echo $? > "$d/rc")
why=
[ "$(cat "$d/rc")" != 0 ] || why=" exit 0 writing to /dev/full"
d=$W/fail.open
mkdir -p "$d"
(cd "$d" && multiz "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf" 0 "$d/no/such/U1.maf" U2.maf > out.maf 2> err; echo $? > rc)
[ "$(cat "$d/rc")" != 0 ] || why="$why exit 0 with an unopenable out1"
result "multiz output errors" "$why"

# a sequence byte of 0x80 or above scores as N; it crashed multiz, and roast hid the crash (findings 14, 2)
d=$W/fail.highbyte
mkdir -p "$d"
cp "$D"/tba/sa.*.sing.maf "$d"
awk 'BEGIN { hi = sprintf("%c", 233) } !done && /^s se/ { n = split($0, f, " "); t = f[n]; i = match(t, /[ACGTacgt]/);
	sub(/[^ ]*$/, substr(t, 1, i-1) hi substr(t, i+1)); done = 1 } { print }' "$D/tba/sa.se.sing.maf" > "$d/sa.se.sing.maf"
roast_run "$d"
why=
cmp -s "$D/tba/sa.se.sing.maf" "$d/sa.se.sing.maf" && why=" no byte was changed"
ran "$d" || why="$why exit $(cat "$d/rc")"
for sp in sd se; do
	grep -aq "^s $sp" "$d/out.maf" || why="$why no $sp rows"
done
result "roast high byte" "$why"

# a crashed child stops roast (finding 2): multiz's first call (the (sd se) merge, after the reference-leaf merge
# that set force) fails; later calls run the real multiz, so only the first failure can be the one reported
d=$W/fail.child
mkdir -p "$d/bin"
cp "$D"/tba/sa.*.sing.maf "$d"
printf '#!/bin/sh\n[ -e "%s/called" ] && exec "%s/multiz" "$@"\n: > "%s/called"\nexit 3\n' "$d" "$B" "$d" > "$d/bin/multiz"
chmod +x "$d/bin/multiz"
(PATH=$d/bin:$PATH; export PATH; roast_run "$d")
why=
[ "$(cat "$d/rc")" != 0 ] || why=" exit 0 after a failed multiz"
result "roast failed child" "$why"

# a headerless MAF whose first 'a' line is long (finding 5)
d=$W/fail.headerless
mkdir -p "$d"
grep -v '^#' "$D/tba/sa.sb.sing.maf" | sed 's/^a score=.*/& extra=padding-padding-padding-padding/' > "$d/in.maf"
(maf_sort "$d/in.maf" sa > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
ran "$d" || why=" exit $(cat "$d/rc")"
result "maf_sort headerless" "$why"

# a 213-character row name comes out whole (finding 9)
d=$W/fail.longname
mkdir -p "$d"
x=$(printf '%0210d' 0)
sed "s/^s sb\([ .]\)/s sb.$x\1/" "$D/tba/sa.sb.sing.maf" > "$d/in.maf"
(maf_sort "$d/in.maf" sa > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
ran "$d" || why=" exit $(cat "$d/rc")"
grep -q "^s sb\.$x " "$d/out.maf" || why="$why name changed"
result "maf_sort long name" "$why"

# a 300-character temporary directory (finding 11)
d=$W/fail.longtmp
mkdir -p "$d"
cp "$D"/tba/sa.*.sing.maf "$d"
lt=$d/$(printf '%0150d' 0)/$(printf '%0150d' 0)
mkdir -p "$lt"
(cd "$d" && roast T="$lt" E=sa "$TREE" sa.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
why=
ran "$d" || why=" exit $(cat "$d/rc")"
grep -v '^#' "$d/out.maf" | cmp -s - "$E/roast.maf" || why="$why out.maf differs"
result "roast long T=" "$why"

# a destination with a space is refused before anything is touched (finding 10)
d=$W/fail.space
mkdir -p "$d/tmp"
cp "$D"/tba/sa.*.sing.maf "$d"
echo keep > "$d/my"
echo keep > "$d/out.maf"
(cd "$d" && roast T="$d/tmp" E=sa "$TREE" sa.*.sing.maf "my out.maf" > /dev/null 2> err; echo $? > rc)
why=
[ "$(cat "$d/rc")" != 0 ] || why=" exit 0"
[ "$(cat "$d/my")" = keep ] && [ "$(cat "$d/out.maf")" = keep ] || why="$why the user's files were changed"
result "roast refuses a destination with a space" "$why"

echo "$pass passed, $fail failed"
[ $fail = 0 ]
