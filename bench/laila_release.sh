#!/usr/bin/env bash
# Release-11.10 timing of scalar (11.9) against the AVX2 build, on laila under hardware counters
# (cannabis-genome scripts/laila_perf.py runs each `run` command under pcount).
#
#   laila_release.sh pack OUT.tar.gz          (local) tar the workloads' inputs, symlinks resolved; KINDS="..." picks
#                                             which (default: deep piece03 cannabis)
#   laila_release.sh build                    (in the job, from the bundle root) build base/ and simd/, unpack $DATA;
#                                             also simd0/ (simd/ built -DYAMA_SIMD=0), simdsh/ (simd/ built
#                                             -DYAMA_SHADOW) and bench/maxrss; EXTRA="dir ..." builds more trees
#                                             of the bundle the same way as base/ and simd/
# BUILD in `run` is base, simd, simd0, simdsh, an EXTRA tree, simdenv0 (simd/'s binaries with YAMA_SIMD=0 in
# the environment: the run-time override a user without AVX2 effectively gets), or simdskip (simd/'s binaries with
# MULTIZ_PREYAMA_SKIP_DASH=1, the experimental pre_yama option). KEEP=1 keeps each run's outputs in its run directory
# (runs/KIND.BUILD.ROUND) instead of deleting them. A simdsh run also prints
#   SHADOW kind round calls=.. rows=.. p1_vec_rows=.. p2_vec_rows=.. unfit_rows=.. bound_fallbacks=..
# summed over every multiz process of the run (a shadow mismatch aborts the run).
#   laila_release.sh run KIND BUILD ROUND     (in the job) one timed run; KIND deep | piece03 | cannabis
#   laila_release.sh compare                  (in the job) every KIND must have one digest across builds and rounds
#   laila_release.sh summarize RESULTS_DIR    (local) per kind and build: median user s, cycles, instructions, IPC
#
# Workloads:
#   deep      the cannabis deep merge, multiz M=1 left.maf right.maf (bench/data/wfi_profile/in)
#   piece03   hg38 100-way piece 03, roast as UCSC's autoMultiz.csh runs it (ucsc100way run03/inputs-v3, tree.nh)
#   cannabis  WF-I piece NW_022060449.1.00: all 216 species, the real tree, run by the multiz_roast UDT's own script
#             (bench/data/wfi_step4/NW_022060449.1), a second workload of a different kind from the human pieces
#   piece04   hg38 100-way piece 04 (ucsc100way run04/inputs-v3), a second human piece
#   cmed      WF-I piece NC_044376.1.02, 216 species, the median-weight piece of the retime study (1.3 GB of inputs)
#   cheavy    WF-I piece NC_044379.1.12, 216 species, the heaviest piece of the retime study (2.2 GB of inputs)
# Each run prints: RESULT kind build round wall_s user_s sys_s digest maxrss_kb (digest: sha256 of the output, '#'
# lines dropped; maxrss_kb: the largest process in the run's tree, from bench/maxrss.c).
set -euo pipefail
B=/home/rcb112/project/claude/bench/data
# bash's own `time` (the gcc image has no /usr/bin/time): wall user sys, written to a file by tm()
TIMEFORMAT='%R %U %S'
tm() { local f=$1; shift; { time "$MAXRSS" "$d/maxrss" "$@"; } 2> "$f"; }
MAXRSS=$PWD/maxrss
# a failed run shows its stderr instead of dying silently under set -e
failrun() { echo "RUN FAILED rc=$1 $kind $build $round" >&2; cat "$d"/err "$d"/*/err 2> /dev/null | tail -20 >&2; exit 1; }
cmd=$1; shift
case $cmd in
pack)
	out=$1 st=$(mktemp -d "$B/tmp/laila_pack.XXXXXX")
	kinds=${KINDS:-deep piece03 cannabis}
	wf() { # KIND CHROM PIECE: a WF-I roast piece's inputs, tree and script
		mkdir -p "$st/$1/inputs"
		cp -L "$B/wfi_step4/$2/pieces/$3"/*.maf "$st/$1/inputs/"
		cp "$B/wfi_step4/$2/roast/$3/tree.nh" "$B/wfi_step4/$2/roast/$3/multiz_roast.sh" "$st/$1/"
	}
	for k in $kinds; do
		case $k in
		deep) mkdir -p "$st/deep"; cp "$B/wfi_profile/in/left.maf" "$B/wfi_profile/in/right.maf" "$st/deep/" ;;
		piece03|piece04) mkdir -p "$st/$k"; cp -L "$B"/ucsc100way/run${k#piece}/inputs-v3/hg38.*.sing.maf "$st/$k/"
			cp "$B/ucsc100way/tree.nh" "$st/$k/tree.nh" ;;
		cannabis) wf cannabis NW_022060449.1 NW_022060449.1.00 ;;
		cmed) wf cmed NC_044376.1 NC_044376.1.02 ;;
		cheavy) wf cheavy NC_044379.1 NC_044379.1.12 ;;
		*) echo "unknown kind $k" >&2; exit 2 ;;
		esac
	done
	tar -C "$st" -cf - $kinds | gzip -1 > "$out"
	rm -rf "$st"
	ls -l "$out"
	;;
build)
	for t in base simd ${EXTRA:-}; do make -s -C "$t" multiz roast maf_project > "build_$t.log" 2>&1 || { tail "build_$t.log"; exit 1; }; done
	rm -rf simd0 && cp -r simd simd0 && make -s -C simd0 clean > /dev/null 2>&1
	make -s -C simd0 multiz roast maf_project CFLAGS="-Wall -Wextra -O2 -DYAMA_SIMD=0" > build_simd0.log 2>&1 || { tail build_simd0.log; exit 1; }
	rm -rf simdsh && cp -r simd simdsh && make -s -C simdsh clean > /dev/null 2>&1
	make -s -C simdsh multiz roast maf_project CFLAGS="-Wall -Wextra -O2 -DYAMA_SHADOW" > build_simdsh.log 2>&1 || { tail build_simdsh.log; exit 1; }
	gcc -O2 -o maxrss simd/bench/maxrss.c
	mkdir -p data && tar -C data -xzf "$DATA"
	ls data/*
	;;
run)
	kind=$1 build=$2 round=$3 bin=$PWD/$build
	if [ "$build" = simdenv0 ]; then bin=$PWD/simd; export YAMA_SIMD=0; fi
	if [ "$build" = simdskip ]; then bin=$PWD/simd; export MULTIZ_PREYAMA_SKIP_DASH=1; fi
	d=$PWD/runs/$kind.$build.$round
	rm -rf "$d"; mkdir -p "$d"
	if [ "$build" = simdsh ]; then export YAMA_SHADOW_LOG=$d/shadow.log; fi
	case $kind in
	deep)
		tm "$d/time" sh -c '"$0" M=1 data/deep/left.maf data/deep/right.maf 0 "$1/U1.maf" "$1/U2.maf" > "$1/out.maf" 2> "$1/err"' \
			"$bin/multiz" "$d" || failrun $?
		dig=$(cat "$d/out.maf" "$d/U1.maf" "$d/U2.maf" | grep -v '^#' | sha256sum | cut -c1-64) ;;
	piece03|piece04)
		mkdir -p "$d/tmp"
		for f in "$PWD"/data/$kind/hg38.*.sing.maf; do ln -s "$f" "$d/"; done
		tree=$(cat data/$kind/tree.nh)
		(cd "$d" && PATH="$bin:$PATH" tm time sh -c 'roast + T="$0/tmp" E=hg38 "$1" hg38.*.sing.maf out.maf > /dev/null 2> err' "$d" "$tree") || failrun $?
		dig=$(grep -v '^#' "$d/out.maf" | sha256sum | cut -c1-64) ;;
	cannabis|cmed|cheavy)
		files=$(ls "$PWD"/data/$kind/inputs/*.maf | paste -sd, -)
		top=$PWD
		(cd "$d" && PATH="$bin:$PATH" tm time sh "$top/data/$kind/multiz_roast.sh" cs10_NCBI_RefSeq_softmasked \
			"$top/data/$kind/tree.nh" "$d/out.maf" "$files") > "$d/stdout" 2> "$d/err" || failrun $?
		dig=$(grep -v '^#' "$d/out.maf" | sha256sum | cut -c1-64) ;;
	*) echo "unknown kind $kind" >&2; exit 2 ;;
	esac
	if [ "$build" = simdsh ]; then
		awk -v k="$kind" -v r="$round" '/^YAMA_SHADOW calls/ { for (i = 2; i <= NF; i++) { split($i, a, "="); v[a[1]] += a[2] } }
			END { printf "SHADOW %s %s calls=%d rows=%d p1_vec_rows=%d p2_vec_rows=%d unfit_rows=%d bound_fallbacks=%d\n",
			      k, r, v["calls"], v["rows"], v["p1_vec_rows"], v["p2_vec_rows"], v["unfit_rows"], v["bound_fallbacks"] }' \
			"$d/shadow.log" | tee -a runs/shadow.txt
	fi
	read -r wall user sys < <(tail -1 "$d/time") # the last line is time's; any earlier ones are the command's stderr
	echo "RESULT $kind $build $round $wall $user $sys $dig $(cat "$d/maxrss" 2> /dev/null || echo NA)" | tee -a runs/results.txt
	[ "${KEEP:-0}" = 1 ] || rm -rf "$d"
	;;
compare)
	rc=0
	for kind in $(awk '{print $2}' runs/results.txt | sort -u); do
		n=$(awk -v k="$kind" '$2 == k {print $8}' runs/results.txt | sort -u | wc -l)
		if [ "$n" = 1 ]; then echo "IDENTICAL $kind $(awk -v k="$kind" '$2 == k' runs/results.txt | wc -l) runs"
		else echo "DIFFER $kind: $n digests"; rc=1; fi
	done
	exit $rc
	;;
summarize)
	# (local) laila_release.sh summarize RESULTS_DIR: per kind and build, medians over rounds of user s, cycles,
	# instructions, IPC and branch-miss rate, joining run_N in counters.tsv to the Nth RESULT line of run.log
	R=$1
	awk 'FNR == NR { if ($1 ~ /^run_/ && $2 != "event") v[$1, $2] = $6; next }
	     /^RESULT/ { n++; r = "run_" n; print $2, $3, $6, v[r, "cycles"], v[r, "instructions"], v[r, "branches"], v[r, "branch-misses"] }' \
		"$R/counters.tsv" "$R/run.log" |
	sort -k1,1 -k2,2 | awk '
	function med(a, m,   t, i, j) { for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (a[j] < a[i]) { t = a[i]; a[i] = a[j]; a[j] = t }
	                                 return m % 2 ? a[(m + 1) / 2] : (a[m / 2] + a[m / 2 + 1]) / 2 }
	function flush() { if (key == "") return
	  printf "%s\t%s\t%d\t%.1f\t%.3e\t%.3e\t%.3f\t%.4f\n", k1, k2, m, med(U, m), med(C, m), med(I, m), med(P, m), med(B, m)
	  split("", U); split("", C); split("", I); split("", P); split("", B); m = 0 }
	BEGIN { print "kind\tbuild\truns\tuser_s\tcycles\tinstructions\tIPC\tbranch_miss_rate" }
	{ if ($1 SUBSEP $2 != key) { flush(); key = $1 SUBSEP $2; k1 = $1; k2 = $2 }
	  m++; U[m] = $3; C[m] = $4; I[m] = $5; P[m] = $5 / $4; B[m] = $7 / $6 }
	END { flush() }'
	;;
*) echo "usage: see the header" >&2; exit 2 ;;
esac
