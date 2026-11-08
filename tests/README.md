# Tests

```bash
make && make check
```

`make check` runs `tests/run.sh`, which needs only `sh`, `awk`, `grep`,
`cmp`, `cp` and `mktemp`. It prints one `PASS` or `FAIL` line per case and a
tally, and exits non-zero if any case fails. It takes about a second.
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

## Expected outputs

The files in `expected/` were made by the 11.6 release. After a change that
is meant to alter the output, regenerate them and review the diff:

```bash
tests/run.sh -g
git diff --stat tests/expected
```

`-g` takes the same optional directory, so expected outputs can come from
another build.

## Inputs

`gen/regen.sh BUILD` rebuilds every input from scratch, then the expected
outputs with that build. It needs python3, and lastz on `PATH` for `all_bz`.
`gen/make_blocks.py` writes the synthetic pairs; `gen/trim.awk` cuts
`k8l8s2` and `k3l40s1` down to two blocks (and file 2 of `k3l40s1` to the
reference and 40 rows, without all-gap columns); `gen/make_species.py`
writes the five species. Run with the 11.6 release and lastz 1.04.52, it
reproduces the committed files byte for byte.

## Not here

The large-data checks need data that is not in the repository and are run
separately: UCSC's hg38 multiz5way rebuilt byte for byte, a whole piece of a
217-genome roast, and pieces of UCSC's hg38 100-way.
