# multiz / TBA

DNA multiple sequence aligner — the Miller Lab release from Penn State, packaged
to build with modern toolchains.

**TBA** (Threaded Blockset Aligner) computes a *threaded blockset*: a
generalisation of a multiple alignment in which every position of every input
sequence appears exactly once, and which can be projected onto any species as a
reference. Different projections of the same blockset are guaranteed consistent
with one another — a property reference-based aligners cannot offer.
**MULTIZ** is the dynamic-programming aligner underneath, and is useful on its
own for fragmented or rearranged genomes.

This software builds the multi-way alignment tracks on the UCSC Genome Browser,
from which their conservation tracks are derived.

## Documentation

| | |
|---|---|
| [**docs/pipelines.md**](docs/pipelines.md) | Concepts, a worked walkthrough, `tba` vs `roast`, genome-scale recipes, troubleshooting. **Start here.** |
| [**docs/programs.md**](docs/programs.md) | Every program, every option, verified against the source. |
| [**docs/formats.md**](docs/formats.md) | Sequence headers, guide trees, the lastz specs file, MAF and filename conventions. |
| [docs/historical/](docs/historical/) | The original upstream `README` and `README2`, unmodified. |
| [tba.pdf](tba.pdf), [tba_howto.pdf](tba_howto.pdf) | The 2004 paper and the 2005 practical guide, as shipped upstream. |

## Building

```bash
make
```

That is the whole build — no configure step, no external libraries. It produces
18 executables in the source directory, and works with current compilers
(gcc 15, whose default is C23).

`make check` then runs the regression tests in a second or so; see
[tests/README.md](tests/README.md).

The Makefile hardcodes `CC = gcc` (`Makefile:1`), so on a system without `gcc`
override it: `make CC=cc`. `make install` additionally needs `install` and
`arch`, and `clean` needs `date`.

To install:

```bash
make && make install INSTALLDIR=/usr/local/bin
```

The `install` target has no prerequisites (`Makefile:77`), so it will not build
for you — run `make` first, as above. `INSTALLDIR` defaults to
`/depot/apps/$(ARCH)/bin`, a Miller Lab convention.

The Makefile compiles at `-O2` with `-Wall -Wextra`. Until 11.3 it used
`-O0`, which made a roast run about twice as slow. For a debug build, override
it: `make CFLAGS="-Wall -Wextra -O0 -ggdb3"`.

## Dependencies

At runtime:

- **[lastz](https://github.com/lastz/lastz)** must be on `PATH`. It is invoked
  by `blastzWrapper`, and therefore by `all_bz`.
- **The installation directory must be on `PATH`.** `tba`, `roast`, and
  `all_bz` invoke other programs in the suite by bare name — see the
  [per-driver table](docs/programs.md#external-programs-not-included).
  `all_bz` runs its steps as `/bin/sh` pipelines, and `blastzWrapper` calls
  `grep`. `roast` and `tba` start their programs directly and do their own
  file handling (since 11.8 and 11.9).

Some non-default options require `toast`, `chain`, or `blastz_clean` from the
separate TOAST distribution, which is **not** included here. The default paths
do not. See
[External programs](docs/programs.md#external-programs-not-included).

## Quick start

```bash
# One sequence file per species, named exactly as in the tree, no extension
for s in human chimp mouse rat chicken; do ln -s data/$s.fa $s; done

TREE="(((human chimp) (mouse rat)) chicken)"

all_bz + "$TREE"                            # 1. pairwise alignments
tba    + "$TREE" *.*.sing.maf tba.maf       # 2. multiple alignment
maf_project tba.maf human > human.tba.maf   # 3. project onto a reference
```

The filename conventions are the interface between the steps and are not
optional — [docs/formats.md](docs/formats.md) explains them.

## The programs

**Pipeline drivers** — `all_bz` runs all the pairwise alignments a guide tree
requires; `tba` builds a reference-independent threaded blockset; `roast`
builds a reference-dependent one from *N*−1 rather than *N*(*N*−1)/2 pairwise
alignments. UCSC drives this last one under the name `autoMZ`.

**Aligners** — `multiz` aligns two ref-blocksets sharing a reference;
`multic` does the same where the reference may contain duplicates;
`blastzWrapper` runs `lastz` over multi-entry FASTA files.

**MAF manipulation** — `maf_project` (project onto a reference), `maf_order`
(reorder and prune species), `maf_sort`, `single_cov2` (enforce single
coverage), `mafFind` (extract an interval), `pair2tb`, `get_covered`.

**Format conversion** — `lav2maf`, `maf2lav`, `maf2fasta`,
`get_standard_headers`.

**Validation** — `maf_checkThread`.

Full option reference: [docs/programs.md](docs/programs.md).

## Provenance and changes from upstream

This repository started as a copy of the multiz/TBA tarball posted at the
Miller Lab website on **21 January 2009** (`http://www.bx.psu.edu/miller_lab`),
released here in 2019 as `20190527` (multiz 11.2). UCSC runs that same release
in production, installed as `multiz.2009-01-21_patched`, with one change: a fix
by **Minmei Hou**, multiz's maintainer in the Miller Lab, which she sent to UCSC
in January 2014. Its source was never published; it was reconstructed from
UCSC's binaries, confirmed by UCSC, and is included here from 11.3 on — see
[`multiz`](docs/programs.md#multiz).

### Changes in 11.9.1

- **Two out-of-bounds reads in `multiz`'s band setup (`pre_yama`) are fixed.**
  Both have been present since at least 2008.
  - One read the 4 bytes before an array wherever a column of the second
    alignment had a reference base and gaps in every other row. Under glibc
    those bytes are always 0, so output was deterministic. Under other
    allocators it was not: under jemalloc 11.9 crashed in 9 of 66 deep merges,
    and under valgrind's allocator it wrote a different cannabis piece. The 0
    is now written explicitly.
  - One passed A's rows as 1..K instead of 0..K-1, so each column was read
    from row 1 to the next column's row 0, and the last column past the end
    of A's buffer.
  - valgrind flagged 201 of 215 `multiz` processes on a 216-genome cannabis
    piece under 11.9 (92 of them for the second read). It reports nothing
    under 11.9.1, and 11.9.1 had no crashes in 66 jemalloc runs.
- **Output can differ in a few blocks** where the old result depended on the
  wrong row range. On two hg38 100-way pieces, 2 and 25 blocks changed (52 and
  421 bp of reference). Aligned bases moved by +9 and −110 and total score fell
  by at most 2 × 10⁻⁶. A deep merge, a 216-genome cannabis piece and
  `make check` (gcc and clang) match 11.9 byte for byte.

Program versions: only `multiz` moves, to 11.9.1. `tba` and `multic` link the
same code, so their output can change in the same rare cases.

### Changes in 11.9

- **No shell for `roast`, almost none for `tba`.** The commands they build
  between merges are now run in-process when they take a simple form:
  programs are started directly, and `mv`, `rm -f`, `touch`, `cat`, `echo` and
  the `grep -v` filters are done without `/bin/sh`. Anything else still goes
  to the shell. On a five-species `tba` run, programs started fall from 203 to
  63 and shells from 99 to 1; `roast` starts no shell.
- **Faster MAF input and output.** `s` lines are written without `fprintf`, and
  lines are read with `getline`. A block read out of file order (as
  `maf_project` reads `roast`'s merged files) is read as exactly its bytes
  instead of refilling a 1 MiB buffer.
- **Failures are reported.** `blastzWrapper` and `lav2maf` stop with an error
  when a lav record is cut short, instead of reusing the previous line.
  `blastzWrapper` reports `lastz`'s failure, so `all_bz` stops. `all_bz`'s
  clean step (TOAST path only) no longer runs in `-` (dry-run) mode, and its
  failure stops the run. `lav2maf` refuses a sequence name too long for its
  buffer instead of overflowing it.
- **Tests and tooling.** `make check` also runs `roast` on an excerpt of UCSC's
  hg38 alignments (24 cases), every test input can be regenerated
  (`tests/gen`), GitHub Actions runs `make check` with gcc and clang, and a
  `.gitignore` lists the built programs.

Output on valid input is unchanged. On the two hg38 100-way pieces, time falls
by 7% against 11.8 (486 s to 453 s, and 496 s to 462 s), and system time from
about 37 s to 25 s. `make check`, the 64-bit shift checks, the deep merge, a
whole 217-genome piece, all seven pieces of an hg38 100-way chr21 rebuild,
UCSC's hg38 multiz5way and a five-species `tba` run all match 11.8 byte for
byte.

Program versions: every program that prints a version moves on — `multiz` 11.9,
`multic` 12.8, `tba` 19, `roast` 8, `maf_project` 18, `lav2maf` 18, `maf2lav`
17, `single_cov2` 17, `pair2tb` 16, `blastzWrapper` 15, `maf_order` 15,
`maf2fasta` 8, `mafFind` 7 and `all_bz` 16.

### Changes in 11.8

- **`roast` without a shell.** It used to run every step between merges —
  `cp`, `mv`, `rm`, `grep` and each `maf_project` and `multiz` — through
  `/bin/sh`, about 2,200 program starts and 1,100 shells on one hg38 100-way
  piece. It now does the file steps itself and starts `maf_project` and
  `multiz` directly, about 300 starts. A leaf's pairwise file is projected
  where it is instead of being copied first.
- **Faster parsing.** MAF `s` lines are read with a dedicated parser instead of
  `sscanf`.
- **`roast` no longer drops rows whose names contain `maf` or `eof`.** The
  steps that remove MAF comment lines were `grep -v maf` and `grep -v eof`,
  which also removed any alignment row whose sequence name contains either
  string, leaving blocks short of a species without an error. Only `#` lines
  are removed now. No name in an hg38 100-way chr21 rebuild or in a
  217-genome *Cannabis* piece contains either string, and their output is
  unchanged.
- **`maf2lav` matches the first species on the minus strand.** It compared the
  start coordinate with `'-'` instead of the strand, so it dropped those
  alignments. `maf2lav` is not used by `roast`, `multiz` or `tba`.
- **`make check`.** Byte-identity tests against committed 11.6 output, and
  coordinate-shift tests; see [tests/README.md](tests/README.md).

Apart from those two fixes, output is unchanged. On the same two hg38 100-way
pieces, time falls by a further 10% and 9% against 11.7 (577 s to 518 s, and
562 s to 509 s), and system time from about 100 s to 40 s; the deep merge,
which `yama` dominates, is unchanged (7.08 s). `make check`, the 64-bit shift
checks, the deep merge, a whole 217-genome piece, all seven pieces of an hg38
100-way chr21 rebuild, UCSC's hg38 multiz5way and a five-species `tba` run all
match 11.7 byte for byte.

Program versions: every program that reads MAF files moves on — `multiz` 11.8,
`multic` 12.7, `tba` 18, `roast` 7, `maf_project` 17, `lav2maf` 17, `maf2lav`
16, `single_cov2` 16, `pair2tb` 15, `blastzWrapper` 14, `maf_order` 14,
`maf2fasta` 7 and `mafFind` 6.

### Changes in 11.7

- **Cheaper dynamic-programming cells.** `yama` merges symbols whose
  substitution scores are identical (upper and lower case, and every
  non-ACGT letter under the shipped matrices), so each cell's substitution
  sum has at most six terms and no branch. The band-edge tests move out of the
  gap sums, and the best of the three node scores is chosen without branches.
- **Larger stream buffers.** MAF inputs and outputs get 1 MiB buffers
  instead of the default (often 4 KiB), which halves the system time of a
  `roast` run. A block read in file order no longer seeks.

Output is unchanged. On one deep merge from a 217-genome roast, time goes from
10.74 s to 7.12 s against 11.6 (medians of three interleaved runs per release
on one host). On two 5-9 Mb pieces of an hg38 100-way chr21 rebuild (99
pairwise inputs), time falls by 33% and 36% (860 s to 577 s, and 880 s to 562
s), and system time by about a third. The 64-bit shift checks, the deep merge,
synthetic merges, a whole 217-genome piece, all seven pieces of an hg38 100-way
chr21 rebuild, UCSC's hg38 multiz5way and a five-species `tba` run all match
11.6 byte for byte.

Program versions: every program that reads MAF files or aligns with `yama`
moves on — `multiz` 11.7, `multic` 12.6, `tba` 17, `roast` 6, `maf_project` 16,
`lav2maf` 16, `maf2lav` 15, `single_cov2` 15, `pair2tb` 14, `blastzWrapper` 13,
`maf_order` 13, `maf2fasta` 6 and `mafFind` 5.

### Changes in 11.6

- **Sequences longer than 2^31 − 1 bases.** Coordinates and sequence lengths
  were 32-bit `int`. In 11.5 a start or length between 2^31 and 2^32 stopped
  with `Size <= 0`, and one past 2^32 wrapped to a wrong coordinate with exit
  status 0. Every program now holds them in 64 bits, including the FASTA and
  lav readers and the sort comparators.
- **Totals and products that overflowed.** `single_cov2`'s count of bases lost
  could pass 2^31 on whole-genome input, and `tba` multiplied two spans as
  `int` when comparing chain connections, which overflowed past about 46,000
  bases each. Both are 64-bit now.
- **Limits that stay 32-bit say so.** A `.nib` file stores a 32-bit length, and
  a block's width is still an `int`. Reading a `.nib` range past that, or a
  `pair2tb` unaligned stretch wider than that, stops with an error instead of
  truncating.
- **No undefined shifts.** `MIN_INT` and `MAX_INT` were built by shifting 1
  into the sign bit. They are now `INT_MIN` and `INT_MAX`.

Output is otherwise unchanged. A test harness shifts every coordinate of real
and synthetic inputs by up to 4 × 10^12 and checks that `roast` and `multiz`
give the same alignment, shifted. 11.6 passes all 165 checks and 11.5 fails 96.
A deep merge, a whole 217-genome piece and UCSC's hg38 multiz5way are byte for
byte the same as 11.5.

Program versions: every program except `all_bz` reads coordinates and moves
on — `multiz` 11.6, `multic` 12.5, `tba` 16, `roast` 5, `maf_project` 15,
`lav2maf` 15, `maf2lav` 14, `single_cov2` 14, `pair2tb` 13, `blastzWrapper` 12,
`maf_order` 12, `maf2fasta` 5 and `mafFind` 4. `get_covered`, `maf_sort`,
`maf_checkThread` and `get_standard_headers` print no version.

### Changes in 11.5

- **Faster block scoring.** `mafScoreRange`, which scores each merged block,
  summed over every pair of rows in every column. It now builds each column's
  sum from symbol counts and four gap classes, as `yama` has since 11.4, and
  the totals are exactly the same. A substitution or gap table that is not
  symmetric under swapping two rows falls back to the pair loop.
- **Less scoring that is never read.** A kept block tail is scored only if it
  is written whole, and a block slice is no longer scored twice.

Output is unchanged. On one deep merge from a 217-genome roast (cs10 and 216
assemblies), time goes from 27.1 s to 11.0 s against 11.4 (medians of three
interleaved runs per release on one host). Every score in that merge and in a
whole 4.4 Mb piece was also checked against the old pair loop. UCSC's hg38
multiz5way still matches byte for byte.

Program versions: every program that calls `mafScoreRange` moves on —
`multiz` 11.5, `multic` 12.4, `tba` 15, `maf_project` 14, `lav2maf` 14,
`maf2lav` 13, `single_cov2` 13 and `mafFind` 3. The others are unchanged.

### Changes in 11.4

- **Faster alignment of deep blocks.** The dynamic-programming step summed its
  gap and substitution scores over every pair of rows in every cell; it now
  builds each sum from per-column counts, which gives exactly the same totals
  (#28). Output is unchanged. On a 216-genome roast, one deep merge goes from
  981 s to 29 s and a whole 4.4 Mb piece from 3.9 h to 0.80 h (medians of
  three interleaved runs per release on one host); UCSC's hg38
  multiz5way still matches byte for byte.

Program versions: `multiz` 11.4, `multic` 12.3 and `tba` 14 (they share the
dynamic-programming code); the others are unchanged from 11.3.

### Changes in 11.3

Alignment output changes only where 11.2 was wrong:

- **Recovered segments.** `multiz` no longer drops the pieces of an overlap it
  cannot align; they go to the unused-block files, which roast re-merges
  (Minmei Hou's 2014 fix, #14). With it, rebuilding UCSC's hg38 multiz5way with
  UCSC's recipe matches their published alignment byte for byte.
- **Long blocks with many rows** are aligned correctly. The dynamic-programming
  scores overflowed 32 bits, which silently misaligned the rest of the block or
  stopped with `Error generating edit script` (#19).
- **Long shared insertions** no longer stop multiz with `ckalloc: request for
  negative space` (#15).
- **MAF blocks run together** without a blank line are read as separate blocks,
  instead of being merged or failing with `inconsistent row size` (#17).

Other changes:

- **Memory.** `multiz`, and `maf_project` as roast calls it, stream their input
  instead of reading it whole: a roast run's peak memory drops from hundreds of
  MiB to a few MiB, with identical output (#24).
- **Speed.** Built at `-O2` by default, about twice as fast as the old `-O0`
  default (#22).
- **MAF files without a header line**, such as `lastz --format=maf-` output, are
  accepted (#21).
- **Temporary files** go to `$TMPDIR` when it is set, and `tba` accepts `T=`
  like `roast` (#18).
- **Robustness.** `roast` and `tba` no longer overflow fixed buffers on long
  command lines (#16), and `formNameList` allocates enough memory (#23).
- **Builds** with current compilers: C23 `bool`, `-fno-common`, strict ISO C
  (#20).

Program versions: `multiz` 11.3, `multic` 12.2, `maf_project` 13, `roast` 4,
`tba` 13; `single_cov2`, `pair2tb` and `maf2lav` 12, `maf2fasta` 4, `mafFind` 2,
`maf_order` 11. Programs whose behaviour did not change keep their numbers.

### Changes in the 2019 release

- **`lastz` replaces `blastz`.** `blastzWrapper` now executes `lastz`
  (`blastzWrapper.c:15`). Blastz has been deprecated for well over a decade and
  lastz is a drop-in replacement. The program keeps its original name so that
  existing scripts and the rest of the suite continue to work. Lastz accepts the
  blastz option letters this suite passes (`Y=`, `H=`, `Q=`), which is what
  makes the substitution work. UCSC appears to have done the same — their TBA
  build notes invoke a `lastzWrapper` — though they publish no source for it.
- Syntax corrections to silence compiler warnings; removal of obsolete `rcsid`
  tags. No behavioural change.
- Addition of the MIT license.

The original motivation for the repository was to provide an equivalent of that
tarball that modern package managers can install. Since 11.3 it is also where
multiz is maintained.

The historical `README` and `README2` are preserved verbatim in
[docs/historical/](docs/historical/). They are a changelog rather than a
reference and have drifted from the code in several places — the current
reference is [docs/programs.md](docs/programs.md), which notes the
discrepancies where they matter.

## Citation

> Blanchette M, Kent WJ, Riemer C, Elnitski L, Smit AFA, Roskin KM, Baertsch R,
> Rosenbloom K, Clawson H, Green ED, Haussler D, Miller W.
> Aligning multiple genomic sequences with the threaded blockset aligner.
> *Genome Research* 14:708–715, 2004. doi:[10.1101/gr.1933104](https://doi.org/10.1101/gr.1933104)

## License

MIT — see [LICENSE](LICENSE). This covers the source code.

The two PDFs in the repository (`tba.pdf`, `tba_howto.pdf`) arrived with the
upstream tarball and are **not** covered by it: `tba.pdf` is the 2004 *Genome
Research* paper, whose copyright rests with its publisher and authors.
