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
for bad in yes 2 0 10; do
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
# multic's pre_yama honours the option too, so its output must record it: data/skipdash, where it changes a block
d=$W/skipdash.multic
mkdir -p "$d"
(cd "$d" && multic "$D/skipdash/hg38.zonAlb1.sing.maf" "$D/skipdash/hg38.taeGut2.sing.maf" 0 > off.maf 2> off.err;
	MULTIZ_PREYAMA_SKIP_DASH=1 multic "$D/skipdash/hg38.zonAlb1.sing.maf" "$D/skipdash/hg38.taeGut2.sing.maf" 0 \
	> on.maf 2> on.err; echo $? > rc)
why=
[ "$(cat "$d/rc")" = 0 ] || why=" exit $(cat "$d/rc")"
grep -q '^# MULTIZ_PREYAMA_SKIP_DASH=1' "$d/on.maf" || why="$why no header line"
! grep -q 'MULTIZ_PREYAMA_SKIP_DASH' "$d/off.maf" || why="$why header line in the default output"
grep -q '^multic[^:]*: MULTIZ_PREYAMA_SKIP_DASH=1' "$d/on.err" || why="$why no notice under multic's name"
result "multic skipdash header" "$why"

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

# full DIR PROGRAM ARGS...: run PROGRAM with stdout on /dev/full; it must fail and say why (a missing /dev/full,
# as on macOS and the BSDs, would make the shell fail instead, so those systems skip it)
full() {
	fd=$1; shift
	mkdir -p "$fd"
	("$@" > /dev/full 2> "$fd/err"; echo $? > "$fd/rc")
	[ "$(cat "$fd/rc")" != 0 ] && grep -q 'error writing' "$fd/err"
}
# unopenable DIR PROGRAM A B: out1, then out2, in a directory that does not exist
unopenable() {
	ud=$1; shift
	mkdir -p "$ud"
	(cd "$ud" && "$@" 0 "$ud/no/such/U1.maf" U2.maf > out.maf 2> err; echo $? > rc1)
	(cd "$ud" && "$@" 0 U1.maf "$ud/no/such/U2.maf" > out.maf 2> err; echo $? > rc2)
	[ "$(cat "$ud/rc1")" != 0 ] && [ "$(cat "$ud/rc2")" != 0 ]
}

# a failed write or an unopenable out1/out2 is an error, not exit 0 (finding 4); multiz and multic
for prog in multiz multic; do
	why=
	if [ -c /dev/full ]; then
		full "$W/fail.full.$prog" $prog "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf" 0 || why=" /dev/full not reported"
	else
		echo "NOTE no /dev/full here: $prog's write-error case is skipped"
	fi
	unopenable "$W/fail.open.$prog" $prog "$D/synth/k2l2s1.1.maf" "$D/synth/k2l2s1.2.maf" || why="$why exit 0 with an unopenable out1 or out2"
	result "$prog output errors" "$why"
done

# a sequence byte of 0x80 or above scores as N; it crashed multiz, and roast hid the crash (findings 14, 2)
d=$W/fail.highbyte
mkdir -p "$d"
cp "$D"/tba/sa.*.sing.maf "$d"
# the byte is made by printf, not awk: gawk and BWK awk write U+00E9 as two UTF-8 bytes in a UTF-8 locale
hi=$(printf '\351')
LC_ALL=C awk -v hi="$hi" '!done && /^s se/ { n = split($0, f, " "); t = f[n]; i = match(t, /[ACGTacgt]/);
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

# 6. Found by the 2026-11-16 code-review-and-quality pass; each case crashed, hung or overran a buffer in 11.9.1
# (and in 11.10 before the fix). None changes the output on valid input.

# an 'a' line of 1000 bytes or more overran a stack buffer in every MAF reader; mafWrite itself writes such
# lines for a block with ~125 paralog rows
d=$W/rev.longa
mkdir -p "$d"
pad=$(printf '%01100d' 0)
sed "s/^a score=\(.*\)/a score=\1 note=$pad/" "$D/tba/sa.sb.sing.maf" > "$d/in.maf"
(maf_sort "$d/in.maf" sa > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
ran "$d" || why=" exit $(cat "$d/rc")"
[ "$(grep -c '^a' "$d/out.maf")" = "$(grep -c '^a' "$D/tba/sa.sb.sing.maf")" ] || why="$why blocks lost"
result "maf_sort 1100-byte a line" "$why"

# an amplifier= or copy= naming a row the block lacks walked off the list (a crash, not a message)
d=$W/rev.amplifier
mkdir -p "$d"
awk '/^a/ && !done { print $0 " amplifier=5"; done = 1; next } { print }' "$D/tba/sa.sb.sing.maf" > "$d/in.maf"
(maf_sort "$d/in.maf" sa > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
[ "$(cat "$d/rc")" = 1 ] || why=" exit $(cat "$d/rc"), not 1"
grep -q 'is not in the block' "$d/err" || why="$why no message"
result "maf_sort amplifier= past the last row" "$why"

# i, e and q lines (UCSC's) are skipped: they leaked memory per line; the output must equal the file without them
d=$W/rev.ieq
mkdir -p "$d"
awk '{ print } /^s sb/ { q = $7; gsub(/[^-]/, "9", q); print "i sb C 0 C 0"; print "e sc 0 10 + 100 I"; print "q sb " q }' \
	"$D/tba/sa.sb.sing.maf" > "$d/in.maf"
(maf_sort "$d/in.maf" sa > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
maf_sort "$D/tba/sa.sb.sing.maf" sa > "$d/want.maf" 2> /dev/null
why=
ran "$d" || why=" exit $(cat "$d/rc")"
grep -q '^i sb' "$d/in.maf" || why="$why no i lines in the input"
grep -v '^#' "$d/want.maf" > "$d/want.data"
grep -v '^#' "$d/out.maf" | cmp -s - "$d/want.data" || why="$why output differs"
result "maf_sort i, e and q lines" "$why"

# FASTA header fields and sequence-file paths were copied unbounded into 200-byte name buffers: pair2tb (which
# tba runs at every merge) aborted on a 350-character path or a 211-character name
d=$W/rev.p2
lp=$d/$(printf '%0170d' 0)/$(printf '%0170d' 0)
mkdir -p "$lp"
cp "$D/tba/sa" "$D/tba/sb" "$D/tba/sa.sb.sing.maf" "$lp"
(pair2tb "$lp/sa.sb.sing.maf" "$lp/sa" "$lp/sb" > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
ran "$d" || why=" exit $(cat "$d/rc")"
result "pair2tb long sequence paths" "$why"

# a 700-species tree: roast built an unused list of every merge's names in name_list[10000] and overran it at
# the root merge, after all the work; a dry run ('-') shows it in a moment
d=$W/rev.bigtree
mkdir -p "$d/tmp"
t=$(awk 'BEGIN { t = sprintf("sp%015d", 0); for (i = 1; i < 700; i++) t = sprintf("(%s sp%015d)", t, i); print t }')
(roast - T="$d/tmp" E=sp000000000000000 "$t" x.maf "$d/out.maf" > "$d/plan" 2> "$d/err"; echo $? > "$d/rc")
why=
[ "$(cat "$d/rc")" = 0 ] || why=" exit $(cat "$d/rc")"
grep -q 'sp000000000000699' "$d/plan" || why="$why the plan never reaches the last species"
result "roast 700-species tree (dry run)" "$why"

# single_cov2 (all_bz's default post-processing) dereferenced NULL when trimming two overlapping blocks left
# both shorter than 6 columns
d=$W/rev.sc2
mkdir -p "$d"
printf '##maf version=1\na score=10\ns sa 100 7 + 1000 ACGTACG\ns sb 100 7 + 1000 ACGTACG\n\na score=10\ns sa 103 5 + 1000 TACGT\ns sb 500 5 + 1000 TACGT\n\n' > "$d/in.maf"
(single_cov2 "$d/in.maf" > "$d/out.maf" 2> "$d/err"; echo $? > "$d/rc")
why=
[ "$(cat "$d/rc")" = 0 ] || why=" exit $(cat "$d/rc")"
result "single_cov2 every block trimmed away" "$why"

# a .nib file (one sequence, ACGTACGT) never reached its end: every read loop read it again forever. head stops
# a looping program by SIGPIPE, so the case cannot hang the suite
d=$W/rev.nib
mkdir -p "$d"
printf '\072\075\351\153\010\000\000\000\041\060\041\060' > "$d/x.nib"
n=$(get_standard_headers "$d/x.nib" 2> /dev/null | head -c 100000 | wc -l)
why=
[ "$n" -le 2 ] || why=" $n lines: the sequence was read again and again"
result "get_standard_headers .nib ends" "$why"

# blastzWrapper ignored lastz's exit status on its reverse and per-sequence paths (the latter pipes through
# grep, whose status pclose saw), and ended the output with #:eof as if lastz had worked. A fake lastz prints a
# little lav and exits with $LASTZ_RC; sx has two sequences, sy one (the reverse path) or two (per sequence)
d=$W/rev.bzw
mkdir -p "$d/bin"
printf '#!/bin/sh\nprintf "#:lav\\nd {\\n  \\"lastz %%s %%s\\n}\\n#:eof\\n" "$1" "$2"\nexit ${LASTZ_RC:-0}\n' > "$d/bin/lastz"
chmod +x "$d/bin/lastz"
printf '>c1\nACGTACGTAC\n>c2\nGGGGCCCCAA\n' > "$d/sx"
printf '>d1\nACGTACGTACGT\n' > "$d/sy1"
printf '>d1\nACGTACGTACGT\n>d2\nTTTTGGGG\n' > "$d/sy2"
why=
for y in sy1 sy2; do
	(cd "$d" && PATH=$d/bin:$PATH blastzWrapper sx $y > ok.$y 2> /dev/null; echo $? > rc.ok.$y)
	(cd "$d" && LASTZ_RC=3 PATH=$d/bin:$PATH blastzWrapper sx $y > bad.$y 2> /dev/null; echo $? > rc.bad.$y)
	[ "$(cat "$d/rc.ok.$y")" = 0 ] && grep -q '^#:eof' "$d/ok.$y" || why="$why $y: a working lastz failed"
	[ "$(cat "$d/rc.bad.$y")" != 0 ] || why="$why $y: exit 0 after lastz failed"
	! grep -q '^#:eof' "$d/bad.$y" || why="$why $y: #:eof after lastz failed"
done
result "blastzWrapper reports a failed lastz" "$why"

# each roast works in its own mkdtemp directory under T= (file names were T= plus the pid, so two runs sharing
# a T= from different containers could clobber each other): two at once in one T= must both give the expected
# output, and T= must be left holding only what was there before
d=$W/rev.workdir
mkdir -p "$d/a" "$d/b" "$d/tmp"
cp "$D"/tba/sa.*.sing.maf "$d/a"
cp "$D"/tba/sa.*.sing.maf "$d/b"
echo keep > "$d/tmp/_MZ_planted"
(cd "$d/a" && roast T="$d/tmp" E=sa "$TREE" sa.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc) &
(cd "$d/b" && roast T="$d/tmp" E=sa "$TREE" sa.*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
wait
why=
for r in a b; do
	ran "$d/$r" || why="$why $r: exit $(cat "$d/$r/rc")"
	grep -v '^#' "$d/$r/out.maf" | cmp -s - "$E/roast.maf" || why="$why $r: out.maf differs"
done
[ "$(ls "$d/tmp")" = _MZ_planted ] && [ "$(cat "$d/tmp/_MZ_planted")" = keep ] || why="$why T= not left as it was: $(ls "$d/tmp" | tr '\n' ' ')"
result "roast work directory" "$why"

# the same with the reference as the RIGHT leaf of its merge, so the force=1 that leaks is the other branch's;
# multiz's first and second calls fail in turn (11.9.1 exits 0 with both). The second call follows a merge of
# two non-reference leaves, which set force too.
for n in 1 2; do
	d=$W/rev.child$n
	mkdir -p "$d/bin" "$d/tmp"
	cp "$D"/tba/sa.*.sing.maf "$d"
	printf '#!/bin/sh\nc=$(cat "%s/n" 2> /dev/null || echo 0); c=$((c + 1)); echo $c > "%s/n"\n[ $c = %s ] && exit 3\nexec "%s/multiz" "$@"\n' \
		"$d" "$d" $n "$B" > "$d/bin/multiz"
	chmod +x "$d/bin/multiz"
	(cd "$d" && PATH=$d/bin:$PATH roast T="$d/tmp" E=sa "((sb sa) (sc (sd se)))" sa.*.sing.maf out.maf > /dev/null 2> err;
		echo $? > rc)
	why=
	[ "$(cat "$d/rc")" != 0 ] || why=" exit 0 after multiz call $n failed"
	result "roast failed child, reference on the right, call $n" "$why"
done

# refused before anything runs: T= with a space or ~, a destination or E= starting with '-', and an argument
# with $ (echoed into the header, it reached /bin/sh). A planted file must survive each.
d=$W/rev.words
mkdir -p "$d/tmp" "$d/sp ace" "$d/t~1"
cp "$D"/tba/sa.*.sing.maf "$d"
echo keep > "$d/-x.maf"
why=
for args in "T=$d/sp ace|E=sa|$TREE|out.maf" "T=$d/t~1|E=sa|$TREE|out.maf" "T=$d/tmp|E=sa|$TREE|-x.maf" \
		"T=$d/tmp|E=-sa|$TREE|out.maf" "T=$d/tmp|E=sa|$TREE \$HOME|out.maf"; do
	(cd "$d" && IFS='|' && set -f && set -- $args && roast "$1" "$2" "$3" sa.*.sing.maf "$4" > /dev/null 2> err; echo $? > rc)
	[ "$(cat "$d/rc")" != 0 ] || why="$why accepted [$args]"
	grep -q 'must' "$d/err" || why="$why no reason given for [$args]"
done
[ "$(cat "$d/-x.maf")" = keep ] || why="$why -x.maf changed"
result "roast refuses unsafe T=, E=, destination and arguments" "$why"

# multi-contig FASTAs with >species.contig headers: parseHeader fell back to using the filename as both
# species and contig, emitting every contig as an unaligned single-row block under the bare species name
# and crashing multiz/get_covered with "Bad coordinates" when a later contig was shorter than an earlier one.
# Also verify pair2tb emits every contig (not just the first) when a pairwise MAF has zero alignments.
d=$W/rev.multicontig
mkdir -p "$d/tmp"
seq1=$(awk 'BEGIN { for (i = 0; i < 40; i++) printf "ACGTACGTAA" }')
seq2=$(awk 'BEGIN { for (i = 0; i < 5; i++) printf "GCTAGCTAGG" }')
for sp in spA spB spC spD; do
	printf '>%s.ctg1\n%s\n>%s.ctg2\n%s\n' "$sp" "$seq1" "$sp" "$seq2" > "$d/$sp"
done
for pair in "spA spB" "spA spC" "spA spD" "spB spC" "spB spD" "spC spD"; do
	set -- $pair
	printf '##maf version=1\na score=1000.0\ns %s.ctg1 10 80 + 400 %s\ns %s.ctg1 10 80 + 400 %s\n\na score=500.0\ns %s.ctg2 5 30 + 50 %s\ns %s.ctg2 5 30 + 50 %s\n\n' \
		"$1" "$(printf '%s' "$seq1" | cut -c11-90)" "$2" "$(printf '%s' "$seq1" | cut -c11-90)" \
		"$1" "$(printf '%s' "$seq2" | cut -c6-35)" "$2" "$(printf '%s' "$seq2" | cut -c6-35)" \
		> "$d/$1.$2.sing.maf"
done
(cd "$d" && tba + T="$d/tmp" "(spA (spB (spC spD)))" sp*.sing.maf out.maf > /dev/null 2> err; echo $? > rc)
why=
ran "$d" || why=" exit $(cat "$d/rc"): $(head -n 2 "$d/err" | tr '\n' ' ')"
grep -q '^s spA\.ctg1 ' "$d/out.maf" || why="$why missing spA.ctg1"
grep -q '^s spA\.ctg2 ' "$d/out.maf" || why="$why missing spA.ctg2"
! grep -qE '^s sp[ABCD] ' "$d/out.maf" || why="$why collapsed bare species row in out.maf"
printf '##maf version=1\n' > "$d/empty.maf"
(cd "$d" && pair2tb empty.maf spA spB > empty_out.maf 2> empty_err; echo $? > empty_rc)
[ "$(cat "$d/empty_rc")" = 0 ] || why="$why pair2tb empty rc=$(cat "$d/empty_rc")"
grep -q '^s spA\.ctg2 ' "$d/empty_out.maf" || why="$why pair2tb dropped second contig on empty MAF"
result "tba and pair2tb multi-contig >species.contig FASTAs" "$why"

echo "$pass passed, $fail failed"
[ $fail = 0 ]
