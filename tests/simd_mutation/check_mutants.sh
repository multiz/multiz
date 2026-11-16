#!/usr/bin/env bash
# Which `make check` fuzz cases (tests/yama_fuzz/check.sh) kill which mutants of the AVX2 kernels.
#   tests/simd_mutation/check_mutants.sh WORKDIR [MUTANT...]     (default: every mutant in mutants.py)
# For each mutant: copy this tree (tracked files, working-tree state), apply it with mutants.py (which exits 1 unless
# it applied exactly once), run check.sh against the copy, and print
#   KILLED name: the FAIL labels      or      SURVIVED name
# check.sh compares a -DYAMA_SIMD=0 build with the dispatch build of the SAME tree, so a mutant in a vector kernel
# or a range check shows up as a difference. The unmutated tree runs first as a control and must pass.
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$H/../.." && pwd)
W=${1:?usage: check_mutants.sh WORKDIR [MUTANT...]}; shift
mkdir -p "$W"
muts=("$@"); [ ${#muts[@]} -eq 0 ] && mapfile -t muts < <(python3 "$H/mutants.py" list)
run() { # NAME [MUTANT]
	local d=$W/$1; rm -rf "$d"; mkdir -p "$d"
	(cd "$SRC" && git ls-files -z | xargs -0 tar -cf - --no-recursion) | tar -C "$d" -xf -
	if [ -n "${2-}" ]; then python3 "$H/mutants.py" apply "$2" "$d/mz_yama.c" > /dev/null || { echo "NOT APPLIED $2"; return; }; fi
	sh "$d/tests/yama_fuzz/check.sh" "$d" > "$d/check.log" 2>&1
	local fails; fails=$(grep '^FAIL' "$d/check.log" | sed 's/^FAIL yama fuzz //; s/ seed [0-9]*//' | sort -u | paste -sd, -)
	if [ -n "$fails" ]; then echo "KILLED $1: $fails"; else echo "SURVIVED $1"; fi
	rm -rf "$d"
}
c=$(run control); echo "$c"
case $c in SURVIVED*) ;; *) echo "control failed: the run is void"; exit 2 ;; esac
for m in "${muts[@]}"; do run "$m" "$m"; done
