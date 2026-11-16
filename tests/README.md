# Tests

```bash
make && make check
```

`make check` runs `tests/run.sh`, which needs only `sh`, `awk`, `grep`,
`cmp`, `cp` and `mktemp`, and then the yama fuzz test (below), which needs
python3. `run.sh` prints one `PASS` or `FAIL` line per case and a tally, and
exits non-zero if any case fails. The two take about 10 s.
`tests/run.sh DIR` tests the programs in `DIR` instead of the source
directory.

Every comparison drops lines starting with `#`, where the programs write
their versions and command lines.

## What each test proves

**`multiz NAME fwd|rev`**: `multiz M=1 A B 0 U1.maf U2.maf` on a synthetic
pair, in both orders, gives the expected `out.maf`, `U1.maf` and `U2.maf`.
The pairs in `data/synth` share a reference row. In `kXlYsN`, file 1 has X
rows and file 2 has Y, from seed N:

| Pair | Blocks | Covers |
|---|---|---|
| `k2l2s1` | 8 | two pairwise inputs |
| `k8l8s2` | 2 | equal depth |
| `k3l40s1` | 2 | unequal depth; one piece goes to the unused-block file |

`k8l8s2` is the first two blocks of a larger pair. `k3l40s1` is the first two
blocks of a 3 × 200 pair with file 2 cut to its first 40 rows and all-gap
columns removed, to keep the test data small.

**`roast`** and **`tba`**: both drivers on lastz alignments of five
synthetic 10 kb species on the tree `((sa sb) (sc (sd se)))`, with
substitutions, indels, two inversions (so some blocks are on the `-` strand),
a moved segment, a duplication and a novel insert. `roast E=sa` uses the four
`sa.*.sing.maf` files, `tba` all ten plus the sequence files. lastz is not
needed: the pairwise alignments are committed, made once by `all_bz` with
lastz 1.04.52.

**`... shifted OFFSET`**: every input row moved along its source by OFFSET
(start + OFFSET, source size + 2 × OFFSET) must give the unshifted output,
moved the same way. The offsets, 2147480000, 5000000000 and 4000000000000,
put coordinates past 2^31 and 2^32. This catches coordinates held in 32
bits; 11.5 fails all twelve. It needs no expected files. `shift.awk` does the
moving.

**`roast name seofb|smafb`**: the `roast` test with species `sb` renamed so
that its name contains `eof` or `maf`, then renamed back, must give the
`roast` output (fields compared, since names set the column padding). roast
drops comment lines containing those strings, and before 11.8 it dropped any
row containing them too. `seofb` covers the merge of the reference with one
leaf, where roast appends the leaf's pairwise file through the `eof` filter;
the old filter loses every `seofb` row. `smafb` guards the `maf` filter, which
sees only blocks a merge left unused, and these inputs leave none, so it
passes with the old filter as well.

**`roast ucsc`** and **`roast ucsc name amafMis1`**: `roast E=hg38` on real
alignments, `data/ucsc`: UCSC's hg38 pairwise alignments with allMis1,
galGal4 and taeGut2 (from their 100-way inputs), cut to the blocks whose hg38
row starts in chr21:5,150,000-5,199,999. The first must give the 11.6 output.
Its merges leave unused blocks, so the second, with allMis1 renamed amafMis1,
exercises the `maf` filter: the old filter keeps 9 of the 13 amafMis1 rows and
writes 23 blocks instead of 24.

**`roast skipdash default|empty|zero|on`**, **`multiz skipdash header`** and
**`multiz skipdash VALUE is off`**: the experimental, opt-in
`MULTIZ_PREYAMA_SKIP_DASH` option. `data/skipdash` is UCSC's hg38 pairwise
alignments with taeGut2 and zonAlb1 (from their 100-way inputs), cut to the
blocks whose hg38 row starts in chr21:16,303,000-16,307,999, where the option
changes one block. With the variable unset, empty or `0`, `roast E=hg38` must give
`roast-skipdash.default.maf` (with `0`, and a warning on stderr); with `1`,
`roast-skipdash.maf` and a notice on stderr. The two files must differ.
`multiz` with `1` must write a `#` header line saying the option is on; with
any other value (`yes`, `2`, `0`) it must warn and give the default output.
Until the 11.10 review those values were fatal, after `multiz` had emptied its
output files, and `roast` exited 0 with an empty alignment.

**Section 5, failures**: each case fails on 11.9.1. `multiz output errors`
(writing to `/dev/full`, and an `out1` that cannot be opened, must exit
non-zero); `roast high byte` (a 0xE9 byte in one row: the sd and se rows must
survive; it crashed `multiz`, and `roast` hid it); `roast failed child` (a
`multiz` whose first call fails must stop `roast`); `maf_sort headerless` and
`maf_sort long name` (a long first `a` line without a header, and a
213-character name, overran buffers); `roast long T=` (a 300-character
temporary directory); and `roast refuses a destination with a space` (which
deleted the user's files `my` and `out.maf`). `run.sh` unsets the variable first, so the other
cases always test the default.

## Expected outputs

The files in `expected/` were made by the 11.6 release. After a change that
is meant to alter the output, regenerate them and review the diff:

```bash
tests/run.sh -g
git diff --stat tests/expected
```

`-g` takes the same optional directory, so expected outputs can come from
another build. It writes nothing unless every case it generates from ran
cleanly, so a failed run can never become the expected output.

## Inputs

`gen/regen.sh BUILD` rebuilds every input from scratch, then the expected
outputs with that build. It needs python3, and lastz on `PATH` for `all_bz`.
`gen/make_blocks.py` writes the synthetic pairs; `gen/trim.awk` cuts
`k8l8s2` and `k3l40s1` down to two blocks (and file 2 of `k3l40s1` to the
reference and 40 rows, without all-gap columns); `gen/make_species.py`
writes the five species. Run with the 11.6 release and lastz 1.04.52, it
reproduces the committed files byte for byte. `data/ucsc` is an excerpt of
UCSC's files, not generated, and `regen.sh` leaves it alone.

## yama fuzz (also run by `make check`)

`tests/yama_fuzz/check.sh` builds yama twice from the tree, as a scalar reference (`-DYAMA_SIMD=0`) and as the
normal build that picks the AVX2 kernels at run time, and requires identical output on random cases from
`gen.py`. It covers plain cases, gap opens of 2^28 that push row coefficients past 32 bits (the scalar fallback),
one range check at a time (`YAMA_FUZZ_GOP_IDX`), substitution scores scaled by 2^23 (`YAMA_FUZZ_SS`, the only
cases that reach the per-row score-sum check), and the run-time override `YAMA_SIMD=0`, also seen from inside a
shadow build (no row may reach AVX2). Both drivers must exit 0, and only stdout is compared, so a crash in both
builds fails. It also tests the
value-bound check that sends a call to the scalar path when its scores could leave the range the AVX2 I scan is
proved for: with `-DYAMA_SHADOW` builds it checks that the check never fires at the shipped limits, and that it
fires, with unchanged output, when each limit is lowered to 0 (`-DYAMA_VBOUND=0`, `YAMA_ZBOUND`, `YAMA_EBOUND`).
It takes a few seconds and needs python3 (it reports SKIP without it). `make check` builds it with make's `CC` and
`CFLAGS`. On a CPU without AVX2, or with `YAMA_SIMD=0` set, both builds run scalar, and it says so. `tests/yama_fuzz/run.sh` compares two trees instead (e.g. against a release).

## CPUs without AVX2

`tests/no_avx2.sh` (needs qemu-user and python3, about 20 s) runs the normal build under
`qemu-x86_64 -cpu Westmere`, a CPU model without AVX2. It checks that yama's fuzz output matches a
`-DYAMA_SIMD=0` build, that a `-DYAMA_SHADOW` build counts no vector rows there (and some under `-cpu Haswell`,
as a positive control), and that the whole `tests/run.sh` suite passes with every program run that way.

## Another allocator

`tests/allocator_check.sh BINDIR [N]` runs `multiz` on a deep merge once
under glibc, then N times under jemalloc (plain, and with
`MALLOC_CONF=junk:true`), and reports crashes and any output that differs
from glibc's. It is not in `make check`, because it needs libjemalloc2
(`apt install libjemalloc2`) and real input data. It exists because two
overruns in `pre_yama` read memory outside their arrays. glibc's chunk
headers made the stray reads harmless, but jemalloc (standing in for macOS
malloc, which bioconda's osx builds use) let stock 11.9 segfault.

## Not here

The large-data checks need data that is not in the repository and are run
separately: UCSC's hg38 multiz5way rebuilt byte for byte, a whole piece of a
217-genome roast, and pieces of UCSC's hg38 100-way.
