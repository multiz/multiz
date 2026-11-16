#!/bin/sh
# Compare yama in two multiz trees on random cases:
#   tests/yama_fuzz/run.sh BASE_TREE TEST_TREE WORKDIR [SEEDS] [NCASES] [MAXDIM] [CFLAGS]
# YAMA_FUZZ_GOP in the environment rescales gap opens in both builds (see main.c).
# Builds main.c against each tree's mz_yama.c and prints MATCH or DIFF per seed.
set -eu
H=$(cd "$(dirname "$0")" && pwd)
BT=$(cd "$1" && pwd) TT=$(cd "$2" && pwd) W=$3 SEEDS=${4:-"1 2 3 4 5"} NC=${5:-500} MD=${6:-200} CF=${7:--O2}
mkdir -p "$W"
for t in base:"$BT" test:"$TT"; do
	n=${t%%:*} d=${t#*:}
	gcc $CF -Wno-free-nonheap-object -I"$d" "$H/main.c" "$d/mz_yama.c" "$d/mz_scores.c" "$d/util.c" -o "$W/yama_$n"
done
rc=0
for s in $SEEDS; do
	python3 "$H/gen.py" "$s" "$NC" "$MD" > "$W/in.$s"
	"$W/yama_base" < "$W/in.$s" > "$W/base.$s"
	"$W/yama_test" < "$W/in.$s" > "$W/test.$s"
	if cmp -s "$W/base.$s" "$W/test.$s"; then
		echo "MATCH seed $s $(sha256sum < "$W/base.$s" | cut -c1-16) $(wc -l < "$W/base.$s") lines"
	else
		echo "DIFF seed $s"; rc=1
	fi
done
exit $rc
