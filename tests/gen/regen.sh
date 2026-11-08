#!/bin/sh
# Rebuild every test input, then the expected outputs, with a trusted build:
#   tests/gen/regen.sh BASE_BUILD_DIR
# Needs python3 and lastz (for all_bz) on PATH. BASE_BUILD_DIR holds all 18
# programs, built from the release the expected outputs should come from.
set -eu
G=$(cd "$(dirname "$0")" && pwd)
T=$(dirname "$G")
B=$(cd "$1" && pwd)
D=$T/data
W=$(mktemp -d "${TMPDIR:-/tmp}/multiz-regen.XXXXXX")
trap 'rm -rf "$W"' EXIT

# synthetic pairs: k2l2s1 whole; k8l8s2 its first 2 blocks; k3l40s1 the first 2
# blocks of k3l200s1 with file 2 cut to the reference and 40 more rows
python3 -I "$G/make_blocks.py" "$W/k2l2s1" 2 2 --blocks 8 --width 1500 --seed 1 > /dev/null
python3 -I "$G/make_blocks.py" "$W/k8l8s2" 8 8 --blocks 8 --width 1500 --seed 2 > /dev/null
python3 -I "$G/make_blocks.py" "$W/k3l200s1" 3 200 --blocks 8 --width 1500 --seed 1 > /dev/null
mkdir -p "$D/synth"
cp "$W/k2l2s1.1.maf" "$W/k2l2s1.2.maf" "$D/synth/"
for i in 1 2; do awk -v nb=2 -v nr=1000 -f "$G/trim.awk" "$W/k8l8s2.$i.maf" > "$D/synth/k8l8s2.$i.maf"; done
awk -v nb=2 -v nr=1000 -f "$G/trim.awk" "$W/k3l200s1.1.maf" > "$D/synth/k3l40s1.1.maf"
awk -v nb=2 -v nr=41 -f "$G/trim.awk" "$W/k3l200s1.2.maf" > "$D/synth/k3l40s1.2.maf"

# five species, aligned by all_bz (lastz)
mkdir -p "$W/sp"
python3 -I "$G/make_species.py" "$W/sp"
(cd "$W/sp" && PATH="$B:$PATH" all_bz "((sa sb) (sc (sd se)))" > /dev/null 2>&1)
rm -rf "$D/tba"
mkdir -p "$D/tba"
cp "$W"/sp/sa "$W"/sp/sb "$W"/sp/sc "$W"/sp/sd "$W"/sp/se "$W"/sp/*.*.sing.maf "$D/tba/"
echo "lastz: $(lastz --version 2>&1 | head -1)"

sh "$T/run.sh" -g "$B"
