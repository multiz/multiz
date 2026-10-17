# Program reference

Every option below was read from the source in this repository rather than
copied from the historical READMEs, which have drifted in several places. Where
the two disagree, the discrepancy is called out under **Notes**. Defaults are
given in parentheses, following the convention of the programs' own usage text.

None of these programs accept `--help`. Run one with no arguments to print its
usage message; every program exits non-zero after printing it.

**Pipeline drivers** — [`all_bz`](#all_bz) · [`tba`](#tba) · [`roast`](#roast)

**Aligners** — [`blastzWrapper`](#blastzwrapper) · [`multiz`](#multiz) · [`multic`](#multic)

**MAF manipulation** — [`maf_project`](#maf_project) · [`maf_order`](#maf_order) · [`maf_sort`](#maf_sort) · [`single_cov2`](#single_cov2) · [`mafFind`](#maffind) · [`pair2tb`](#pair2tb) · [`get_covered`](#get_covered)

**Format conversion** — [`lav2maf`](#lav2maf) · [`maf2lav`](#maf2lav) · [`maf2fasta`](#maf2fasta) · [`get_standard_headers`](#get_standard_headers)

**Validation** — [`maf_checkThread`](#maf_checkthread)

---

## Pipeline drivers

### all_bz

Generates and runs all pairwise alignments needed to seed a multiple alignment.

```
all_bz [-+] [b=?] [A=?] [F=reference] [T=annotation-file] [h=?] [q=?] [s=?] [c=?] [D=?] [f=?] species-guide-tree [specs-file]
```

| Option | Meaning |
|---|---|
| `+` | Verbose: echo each command to stdout before running it. |
| `-` | Dry run: echo commands only, do not execute. |
| `b` (2) | `0` post-process only · `1` run the aligner only · `2` both. |
| `A` (1) | Post-processor: `0` toast · `1` `single_cov2` · `2` toast, then chain and single-coverage on the reference. |
| `F` (null) | Reference species. Affects only pairs that *contain* the reference: those get single coverage on the reference alone. Pairs of two non-reference species still get reciprocal single coverage (`all_bz.c:230`). |
| `D` (1) | `0` generate the pairs `roast` needs · `1` generate the pairs `tba` needs. |
| `T` (null) | Annotation file, passed to `toast` as `A=<file>`. |
| `h` (300) | Minimum chaining size — `toast` only. |
| `q` (600) | Minimum cluster size — `toast` only. |
| `s` (0) | `0` full-size toast · `1` reduced-size toast. Absent from the program's own usage text; the default is the initialiser at `all_bz.c:61`. |
| `c` (500) | Passed to `blastz_clean`; alignments closer than this are subject to cleaning. |
| `f` (2) | Percentage used to determine in-paralogs — `toast` only. |

The species pairs are derived from the guide tree. `D=1` (the default) produces
the *N*(*N*−1)/2 cross pairs that `tba` needs — or *N*(*N*+1)/2 including self
pairs when `A=0` (`all_bz.c:394`). `D=0` produces the *N*−1 pairwise files that
`roast` needs and requires `F=`; note it additionally runs *N* **self**
alignments (`all_bz.c:384,392`), so the aligner is invoked 2*N*−1 times even
though only *N*−1 files are post-processed.

For each pair, `all_bz` runs (`all_bz.c:46,48-49`):

```
blastzWrapper sp1 sp2 Y=9000 H=0 <specs-options> \
  | lav2maf /dev/stdin sp1 sp2 \
  | maf_sort /dev/stdin sp1 > sp1.sp2.orig.maf
single_cov2 sp1.sp2.orig.maf [R=ref] > sp1.sp2.sing.maf
```

`Y=9000 H=0` is hardcoded and always applied; options from the specs file are
appended to it.

**Notes**

- `README2` documents this option as `t(0)`; it is now **`s=`**. The internal
  error message still says "argument t", which is a leftover.
- `T=`, `c=`, and `f=` are accepted by this version but absent from `README2`.
  `A=2` is likewise newer than the documented `A=0|1`.
- `A=0` requires `blastz_clean` and `toast`; `A=2` additionally requires
  `chain`. None are [part of this package](#external-programs-not-included).
- Version 15 (`all_bz.c:41`).

### tba

Reference-independent multiple aligner. Builds a threaded blockset by walking
the guide tree from leaves to root, invoking `multiz` (or `multic`) at each
internal node.

```
tba [+-] [R=?] [M=?] [E=?] [P=?] [X=?] species-guide-tree maf-source destination-file
```

| Option | Meaning |
|---|---|
| `+` | Verbose: echo each sub-command before running it. |
| `-` | Dry run: echo commands only. |
| `R` (30) | Dynamic-programming radius, passed to the aligner. |
| `M` (1) | Minimum block length in the output. |
| `E` (null) | With `E=species`, produce a reference-centric alignment with single coverage guaranteed for that species only. Without it, single coverage is guaranteed for every species. |
| `P` (null) | Aligner to use: default `multiz`; `P=multic` selects `multic`. |
| `X` (0) | Which pairwise suffix to consume: `0` `.sing.maf` · `1` `.toast.maf` · `2` `.toast2.maf`. |

`maf-source` is the list of pairwise MAF files, normally given as a glob such as
`*.*.sing.maf`. `destination-file` is required.

**Undocumented alternative:** `tba [opts] -f listfile tree destination` reads the
MAF filenames from `listfile`, one per line, instead of from the command line
(`tba.c:387`). `roast` accepts the same form (`auto_mz.c:237`). This is useful
when the list would overflow the shell's argument limit.

**`tba` requires the sequence files to be present** in the working directory,
named exactly by species name — it passes those names to `pair2tb` verbatim
(`tba.c:129`). Supplying only the MAF files is not enough.

**`tba` shells out** to `multiz`, `multic`, `maf_project`, `pair2tb`, and
`get_covered`, plus `cat`, `echo`, `grep`, `mv`, `rm`, and `touch`. All of these
must be on `PATH`. (`roast` shells out to a smaller set — only `multiz`,
`multic`, `maf_project`, plus `cp`, `echo`, `grep`, `mv`, `rm`.) See
[External programs](#external-programs-not-included).

**Notes**

- If `E=` is set, single coverage is only enforced on the reference, so the
  inputs may contain duplicates relative to other species — use `P=multic`.
- `tba` does **not** accept `T=` (temp directory); only [`roast`](#roast) does.
- Version 12 (`tba.c:5`).

### roast

Reference-dependent multiple aligner. Same progressive tree walk as `tba`, but
every pairwise alignment not involving the reference is discarded, so it needs
only *N*−1 pairwise alignments instead of *N*(*N*−1)/2 — a large reduction in
running time.

```
roast [+-] [R=?] [M=?] [P=?] [T=?] [X=?] [C=?] E=reference-species species-guide-tree maf-source destination
```

| Option | Meaning |
|---|---|
| `+` / `-` | Verbose / dry run. |
| `R` (30) | Dynamic-programming radius. |
| `M` (1) | Minimum block length of output. |
| `P` (multiz) | `multiz` requires single coverage on the reference row; `multic` imposes no such requirement. |
| `T` (/tmp) | Alternate temporary directory. |
| `X` (0) | Pairwise suffix selector, as for `tba`. |
| `C` (50) | Connection threshold, 0–100. **Only valid with `P=multic`.** |
| `E` | Reference species. **Required** — `roast` exits if it is absent. |

**Notes**

- **UCSC drives this program under the name `autoMZ`.** Their multi-way
  alignment recipes call `autoMZ` out of a directory named
  `multiz.2009-01-21_patched` — the same upstream tarball this repository
  forked — and invoke it as `autoMZ + T=$tmp E=$db "tree" *.sing.maf out`,
  an argument set that only `auto_mz.c` accepts (`tba` has no `T=`). The name
  `roast` was introduced when it replaced `autoMZ` upstream on 2008-05-07
  (`docs/historical/README2:2`). Treat UCSC's `autoMZ` recipes as `roast`
  recipes; the correspondence is by interface, as their build is patched and
  was not compared byte for byte.
- `C=` is forwarded to the aligner unconditionally. Because `multiz` only scans
  for `RMLS` (`multiz.c:205`), a leading `C=` ends its option loop immediately
  and the argument count check at `multiz.c:238` then fails — so `roast C=…`
  without `P=multic` dies with `multiz`'s **usage dump**, not a message naming
  `C`. (`multiz.c:223`'s `illegal flag` branch is unreachable.)
- `T=` is accepted by this version but absent from `README2`. UCSC relies on it
  to place scratch files on `/dev/shm`.
- `roast` is stated to guarantee that every input pairwise alignment appears in
  the final blockset, which `tba` cannot guarantee, and that every non-reference
  row is orthologous to the reference (*partial orthology*). All pairs of rows
  are orthologous only when the guide tree is skewed and the reference sits at
  the innermost node. These are the algorithm's designed properties as described
  in Hou (2007), §4.2.4 — they are not verified by anything in this repository.
- Internal version 3 (`auto_mz.c:7`).

---

## Aligners

### blastzWrapper

Runs the pairwise aligner over two multi-entry FASTA files.

```
blastzWrapper seq-file1 seq-file2 [options]
```

All `[options]` are passed straight through to the underlying aligner.

**In this fork the underlying aligner is `lastz`, not `blastz`**
(`blastzWrapper.c:14`). `lastz` must be on `PATH`. The program name is
unchanged so that existing scripts and the rest of the suite keep working.

The wrapper exists because the aligner only handles the first FASTA entry of
its first input. `blastzWrapper` runs it repeatedly so that every entry of one
file is aligned against every entry of the other. Note that it **swaps the two
files** when the first holds more entries than the second, iterating over the
smaller one (`blastzWrapper.c:101-111`); the set of pairs is unaffected, but
"file 1" internally is not necessarily the file you passed first. It also
shells out to `grep` (`blastzWrapper.c:132`).

Each entry must carry a
[formatted header](formats.md#formatted-headers-multiple-contigs-per-species)
when a file holds more than one.

Output is in `lav` format on stdout; pipe it through [`lav2maf`](#lav2maf).

Version 11 (`blastzWrapper.c:15`).

### multiz

Aligns two ref-blocksets that share the same reference, guided by a pairwise
blockset. This is the dynamic-programming core of `tba` and `roast`.

```
multiz [R=?] [M=?] [L=?] [S=?] file1 file2 v [out1 out2] [nohead] [all]
```

| Option | Meaning |
|---|---|
| `R` (30) | Radius: the guiding alignment's path through the DP grid is widened by this many cells on each side. |
| `M` (1) | Minimum output width; `1` outputs all blocks. |
| `L` (20) | Large break width. **Undocumented upstream, and inert** — see below. |
| `S` (2) | Small break width. **Undocumented upstream, and inert** — see below. |
| `v` | Required. `0` neither reference alignment is fixed; `1` the reference alignment in `file1` is fixed. |
| `out1 out2` | Files collecting unused blocks from `file1` and `file2`. Both or neither; default is stdout. **Incomplete in this build** — see below. |
| `nohead` | Suppress the MAF header. |
| `all` | **Inert in this build** — accepted and consumed, but single-row blocks are emitted either way. See the note below. |

`file1` and `file2` are MAF files each topped by the same reference sequence.
**Neither may contain duplicates of the reference** — use [`multic`](#multic) if
they might.

`L=` and `S=` parse and validate, but **nothing reads them**. The only code that
used them — deciding when to leave a block unbroken rather than split it — sits
inside a comment block spanning `multiz.c:99-115`, headed "the unused positions
overlap with the unbroken part". Their defaults live in `multi_util.c:15-16`;
setting them changes no output. They appear in neither the usage message nor any
historical README.

**`out1`/`out2` do not capture everything.** `multiz` passes only `fpw2` into
`pre_yama()` (`multiz.c:149`), so when that routine takes one of its early
returns the corresponding segment of `file1` is written nowhere — and on two of
the three bail-out paths the `file2` segment is dropped as well. Only the
`K==0` path emits anything (`mz_preyama.c:195`). Blocks lost this way appear in
neither the alignment nor the unused-block files.

This is what UCSC's patch to their `multiz.2009-01-21_patched` build addresses.
The patch itself is unpublished, but it has been reconstructed: the
reconstruction threads a second `FILE* fpw1` through `pre_yama` and adds
`print_part_ali_col` calls at the two silent bail-outs, so the unaligned
segments are recovered into both unused streams. If you compare this tree's
`out1`/`out2` output against UCSC's, expect theirs to contain strictly more.

The reconstruction has been verified at the binary level. Rebuilt with the two
compilers UCSC used (gcc 4.4.6-4 and 4.4.7), `multiz`, `roast`, and
`maf_project` come out at UCSC's exact file sizes with byte-identical machine
code and data — `.text`, `.rodata`, `.data`, `.init`, `.fini`, `.plt`,
`.eh_frame`, and the dynamic symbols all match. The only differences are 326
bytes per binary: the random build ID and the compiler-runtime version note.

So for those three programs the unused-block change is not merely correct where
it overlaps — it is the *only* change. Two programs fall outside that evidence:
`tba`, which UCSC does not use in their multiz builds, and `multic`, which they
do not ship at all.

**If you apply that patch here, keep the `K==0` NULL guard.** UCSC's build drops
it, which is safe for them because they ship no `multic` and `roast` defaults to
`P=multiz`, so nothing ever passes a null stream. In this tree `multic` is built
and reachable — via `tba P=multic` or `roast P=multic` — and `multic.c:72` is
the one caller that passes `NULL`, which neither `print_part_ali_col`
(`multi_util.c:620`) nor `mafWrite` (`maf.c:251`) checks for.

This is not a theoretical edge case. Tested on three ordinary ~2 MB buckets of
the hg38 5-way under `P=multic`:

| `multic` build | Result |
|---|---|
| stock 2009-01-21 | runs — the reference output |
| patch **with** the guard | identical alignment lines to stock on all three |
| patch **without** the guard (UCSC's) | segfaults on all three — `fprintf(NULL)`, then `maf_project` fails |

The guarded form is the one to use: identical to UCSC's under `multiz`,
identical to stock under `multic`.

**`all` does nothing in this build.** The flag sets `row2=0` (`multiz.c:229`),
which is already the initial value (`multi_util.c:24`), and every output guard
reads `if (row2==0 || ali->components->next != NULL)`
(`multiz.c:70,76,279,284`). The left disjunct is therefore always true and
single-row blocks are **always** written, with or without `all`. The source
comment at `multiz.c:3` ("make it default not to output") records the intent,
but the code does not implement it. The same applies to `multic`
(`multic.c:319,389,394`). Note that [`maf_order`](#maf_order) is unaffected — it
uses its own local flag and does behave as documented.

**Notes**

- Version 11.2 (`multiz.c:58`). `README2` describes v10.6.
- The published default radius is `R = 10` and this build defaults to `R = 30`.
  Blanchette et al. (2004) report, in the Methods supplement ("Aligning
  Procedures"), that an evaluation by evolutionary simulation found `R = 10`
  works essentially as well as `R = 50` — so the choice is not very sensitive.

### multic

Identical interface to `multiz`, but tolerates duplicates of the reference in
both input files.

```
multic [R=?] [M=?] [C=?] file1 file2 v [out1 out2] [nohead] [all]
```

`R`, `M`, `v`, `out1 out2`, `nohead`, and `all` behave exactly as in
[`multiz`](#multiz). Additionally:

| Option | Meaning |
|---|---|
| `C` (50) | Connection threshold, 0–100. **Inert in this build** — see below. |
| `s` (0) | Alignment category, controlling how rows marked as paralogous are handled (`multic.c:133-171`). Undocumented upstream. |

Use `multic` whenever single coverage has not been enforced on every species —
in particular with `tba E=species`, where only the reference is guaranteed
single-coverage.

**`C=` has no effect in this build.** It is range-checked and assigned to
`CONNECTION_THRESHOLD` (`multic.c:307`), but nothing reads that variable at
runtime. Of its two uses, `align_util.c:510` lies inside a comment block
spanning `align_util.c:406-519`, and `align_util.c:656` is inside
`connectionAgreement2()`, whose sole caller is `pre_yama2()`
(`mz_preyama.c:436`) — and `pre_yama2` is never called by anything
(`mz_preyama.c:387` is its only definition, with a prototype at
`mz_preyama.h:27`). `multic` itself calls `pre_yama()` (`multic.c:72`).

**Notes**

- `multic` does **not** accept `L=` or `S=`.
- Version 12.1 (`multic.c:35`). `README2` describes v11.1.

---

## MAF manipulation

### maf_project

Extracts the blocks naming a given reference and orders them along it. This is
how a threaded blockset becomes a conventional reference-based alignment.

```
maf_project file.maf reference [from to] [filename-for-other-mafs] [species-guide-tree] [nohead]
```

| Argument | Meaning |
|---|---|
| `reference` | Sequence to project onto. It becomes the top row of every output block, and blocks are ordered by its start position. |
| `from to` | Restrict output to this interval of the reference. |
| `filename-for-other-mafs` | Collect blocks not present in the projection. **When this is omitted the projected blocks are "beautified"** — adjacent blocks are fused where possible. |
| `species-guide-tree` | Species not named here are removed from the output. |
| `nohead` | Suppress the MAF header. |

Handles species with more than one contig.

**The optional arguments are positional, not freely combinable.** Parsing is
driven by the argument count (`maf_project.c:574-608`), with these consequences:

- `maf_project f.maf ref 0 1000 other.maf` does **not** do what it looks like.
  With five arguments the last one is taken as a *species tree*
  (`maf_project.c:582`), so `other.maf` is parsed as a species list, no
  other-mafs file is written, and the output is filtered to a species named
  `other.maf` — an empty projection, with no warning.
- `filename-for-other-mafs` and `species-guide-tree` can never be given
  together; that argument count is rejected outright.
- `from to` together with both trailing options is likewise rejected.

In practice, use at most one of the three optional trailing arguments at a
time. `nohead` is stripped first (`maf_project.c:574`) and so may always be
last.

Projecting does not bias the alignment toward the chosen species — it is still
the same multiple alignment, viewed from one sequence. Any two projections of
the same threaded blockset are guaranteed consistent: if position *x* of *X*
aligns to position *y* of *Y* in one projection and to position *z* in another,
then *y* = *z*.

Version 12 (`maf_project.c:29`). `README2` describes v10.

### maf_order

Reorders the rows of every block to a given species order, dropping species not
named.

```
maf_order maf-file species1 species2 ... [nohead] [all]
```

| Argument | Meaning |
|---|---|
| `species1 species2 ...` | Rows appear in this order. Species not listed are excluded. |
| `nohead` | Suppress the MAF header. |
| `all` | Keep single-row blocks; by default they are dropped. |

Useful for pruning species from a finished alignment without re-running `tba`.

### maf_sort

Sorts blocks by the position of a named species' row.

```
maf_sort maf-file species-name [unused-ali-file]
```

`[unused-ali-file]` collects blocks that do not contain the named species.
**If you omit it, those blocks are discarded silently** (`maf_sort.c:50`) — no
warning, no count.

Sorting is not all it does: the named species' row is promoted to the top of
each block, and if that row is on the minus strand the whole block is reverse
complemented (`maf_sort.c:34-43`). Used inside the `all_bz` pipeline to
normalise freshly converted pairwise alignments.

### single_cov2

Removes overlapping regions so that each position of a sequence participates in
at most one alignment — "single coverage".

```
single_cov2 pairwise.maf [R=species] [F=deleted.maf]
```

| Option | Meaning |
|---|---|
| `R` (null) | Enforce single coverage on this species only. By default it is enforced reciprocally on both. |
| `F` (null) | File collecting the removed regions. |

Blocks need not be sorted, but the first rows of all blocks must belong to one
species and the second rows to another — run [`maf_project`](#maf_project) or
[`maf_order`](#maf_order) first if that is not already true.

**Notes**

- **Argument order matters.** The parser inspects the last argument for `F=`,
  then the new last argument for `R=` (`single_cov2.c:180-186`). Write
  `single_cov2 in.maf R=human F=del.maf`. Reversing them causes `F=` to be
  **silently ignored** — no warning, no deleted-regions file.
- The usage text says "if S=species specified"; the flag is actually `R=`.
- The `2` in the name means it operates on MAF rather than `lav`.

### mafFind

Prints the blocks intersecting an interval.

```
mafFind file.maf beg end [species-prefix] [slice]
```

| Argument | Meaning |
|---|---|
| `beg end` | Interval to search. |
| `species-prefix` | Match against *any* row whose `src` starts with this prefix, so `mm3` matches `mm3.chr7`. The first row is also a candidate (`mafFind.c:54-57`), so a prefix matching the reference will match there too. |
| `slice` | Trim the reported blocks so their ends match `beg`–`end` exactly. |

With no `species-prefix`, blocks whose **first** row intersects the interval are
printed.

### pair2tb

Converts a pairwise MAF into a threaded blockset.

```
pair2tb pairwise.maf seq-file1 seq-file2
```

Assumes the input has no overlapping blocks, that the top rows correspond to
`seq-file1` and the second rows to `seq-file2`. Both files may contain multiple
contigs. Called by `tba` at every leaf pair.

### get_covered

```
get_covered file1 file2
```

Prints the portions of the blocks in `file1` that are covered by blocks in
`file2`, on stdout. An internal helper invoked by `tba`; it has no usage message
beyond `arguments: file1 file2` and is rarely run by hand.

---

## Format conversion

### lav2maf

```
lav2maf lav-file seq-file1 seq-file2
```

Converts aligner output in `lav` format to MAF. `seq-file1` and `seq-file2` are
the same sequence files given to the aligner — they are needed to recover the
sequence text and header metadata. Accepts `/dev/stdin`, which is how `all_bz`
uses it.

Version 13 (`lav2maf.c:3`).

### maf2lav

```
maf2lav align.maf seq1 seq2
```

The inverse of `lav2maf`. Produces `lav` output suitable for viewers such as
`laj`.

### maf2fasta

Converts a projected MAF into FASTA or MultiPipMaker text-alignment format,
representing the entire reference sequence but only the aligning portions of
the others.

```
maf2fasta refseq-file maf-file [beg end] [fasta[2]][?] [iupac2n] [refsrc=src]
```

| Argument | Meaning |
|---|---|
| `beg end` | Limit output to this interval of the reference. |
| `fasta` / `fasta2` | Request FASTA output instead of the default MultiPipMaker format. With `fasta2`, rows are wrapped at 50 columns (`COL_WIDTH`). |
| `fasta@` (any suffix char) | Append a character to `fasta`/`fasta2` and it replaces `-` **between** two local alignments, distinguishing between-alignment gaps from within-alignment gaps. |
| `iupac2n` | Map any nucleotide outside `ACGTNacgtn` to `N` or `n`, preserving case. Applies **only to the reference sequence** read from `refseq-file` (`maf2fasta.c:98-99`); rows taken from the MAF are untouched. |
| `refsrc=src` | Set the reference species' component source explicitly. Defaults to the source of the first component of the first block. |

The input **must be projected onto the reference** and must contain no
inversions or overlaps. Run [`maf_project`](#maf_project) first.

### get_standard_headers

```
get_standard_headers seq-file
```

For each entry it prints two lines (`get_standard_headers.c:27-28`):

```
<the entry's existing header> ==>
1-<len>:+:<len>
```

Note the **hyphen** between start and end, and that `end` and `srcSize` are the
same value — both are the sequence length. The output is a report, not a
paste-ready header: you still have to compose
`species:contig:start:strand:srcSize` yourself.

`docs/historical/README:19-23` describes this as emitting `1:end:+:srcSize`
with `srcSize` "one larger than end". **That is not what this version does** —
there is no off-by-one, and the separator is not a colon.

---

## Validation

### maf_checkThread

```
maf_checkThread projected-maf-file
```

Checks the threading condition on the **first row** of each block: that the
reference's rows tile its sequence without overlap, in increasing order. Prints
one line per violation and a `Total Errors: N` summary.

- The input must already be projected onto the reference.
- **Only overlap is detected, not gaps.** The sole test is
  `plusStart < lastEnd + 1` (`maf_checkThread.c:26`), so blocks that leave a
  hole in the reference pass silently. Strict threading requires no gaps; this
  program does not check for them.
- Violations are printed **without a trailing newline**
  (`maf_checkThread.c:27`), so several run together on one line.
- Only the first row is checked. To confirm a TBA alignment is fully threaded
  you must project onto each species in turn and check each projection.
- `docs/historical/README:253` says a sequence not starting at 0 yields one
  spurious error. **It does not**: `lastEnd` starts at `-1`
  (`maf_checkThread.c:21`), so the first block's test is `plusStart < 0`, which
  a non-negative MAF coordinate never satisfies. Such a projection reports zero
  errors.

**Note:** `README2:96` documents an `[overlap]` argument. **It does not exist in
this version** — `maf_checkThread.c` reads `argv[1]` only and ignores anything
further.

---

## External programs not included

Several options invoke programs that are **not part of this package** and are
not built by the Makefile. Selecting them without a separate installation fails
at the point the sub-command runs.

| Program | Required by | Source |
|---|---|---|
| `lastz` | `blastzWrapper`, and therefore `all_bz` | [lastz](https://github.com/lastz/lastz) |
| `toast` | `all_bz A=0`, `A=2`; `tba`/`roast` `X=1` | TOAST distribution |
| `chain` | `all_bz A=2`; `tba`/`roast` `X=2` | TOAST distribution |
| `blastz_clean` | `all_bz A=0`, `A=2` | TOAST distribution |
| `laj` | **Never invoked.** The only call site is inside `#ifdef MAF2LAV` (`all_bz.c:419`), and `MAF2LAV` is defined nowhere in the tree or the Makefile. Listed here only because the string appears in the source. | — |

The default paths — `all_bz A=1` and `tba`/`roast X=0` — need only `lastz`.

In addition, three drivers invoke other programs **from this package** by bare
name, so the installation directory must be on `PATH` — not merely the current
directory:

| Driver | Invokes from this package |
|---|---|
| `tba` | `multiz`, `multic`, `maf_project`, `pair2tb`, `get_covered` |
| `roast` | `multiz`, `multic`, `maf_project` |
| `all_bz` | `blastzWrapper`, `lav2maf`, `maf_sort`, `single_cov2` |

Between them they also use the standard utilities `cat`, `cp`, `echo`, `grep`,
`mv`, `rm`, and `touch`.

---

## Source files that are not built

Two C files in the repository are not referenced by the Makefile and produce no
binary:

- `maf_project_simple.c` — declares `VERSION 20` and carries the *same* usage
  string and argument parser as `maf_project` (`maf_project_simple.c:19-20`);
  what it lacks is the block-fusing "beautify" pass. Removed from the build on
  2004-10-25 (`docs/historical/README:30`).
- `dna_nib.c` — unreferenced.

They are retained for provenance with the upstream tarball.
