# Running an alignment

- [Concepts](#concepts)
- [The three steps](#the-three-steps)
- [Walkthrough: a small TBA alignment](#walkthrough-a-small-tba-alignment)
- [tba or roast?](#tba-or-roast)
- [Driving the pipeline by hand](#driving-the-pipeline-by-hand)
- [Genome-scale alignment: how UCSC does it](#genome-scale-alignment-how-ucsc-does-it)
- [Troubleshooting](#troubleshooting)

## Concepts

The vocabulary below comes from Blanchette et al. (2004) and is used precisely
throughout the code and this documentation.

**Block** — a rectangular array of nucleotides and dashes, where deleting the
dashes from any row yields a run of consecutive positions in one of the original
sequences or its reverse complement. No column is all dashes. **A block may
have a single row**, representing sequence that participates in the blockset but
aligns to nothing.

**Blockset** — a set of blocks.

**Threading** — a sequence *S* *threads* a blockset if every position of *S*
appears in exactly one block. Not "is present somewhere", but *exactly once*:
the blocks tile *S* without gaps or overlaps.

**Threaded blockset** — a blockset threaded by *every* original sequence. This
is what `tba` produces, and it is the generalisation of a multiple alignment
that gives the suite its name.

**Ref-blockset** — a blockset where *every* block has a designated row, all
from the same sequence (the *reference*), and each reference position appears in
exactly one block. Note that a threaded blockset is **not** itself a
ref-blockset for each of its sequences: a match need not involve every sequence,
so some blocks lack any row from a given species. That is exactly why a
ref-blockset has to be *generated*.

**Projection** — generating a ref-blockset from a threaded blockset by keeping
the blocks that contain the chosen species, **discarding those that do not**,
and ordering the rest along it. This is what
[`maf_project`](programs.md#maf_project) does.

The property that makes all of this worth the trouble: **projections never
contradict each other.** In the paper's words, if position *x* of *X* aligns to
position *y* of *Y* in one projection and to position *z* of *Y* in another,
then *y* = *z*. Note what this does and does not promise: *x* may align to
nothing at all in a third projection, because the block holding it was dropped
for lacking that species. What is ruled out is two projections making
*conflicting* claims. Reference-based aligners cannot offer even that — align
everything to human and everything to mouse separately, and the two answers can
genuinely disagree.

### What TBA does not do

The threaded-blockset *concept* accommodates inversions and duplications, but
**the `tba` program does not produce them**. It is restricted to matches
occurring in the same order and orientation in all sequences. Concretely, no
output block contains a row on the reverse strand, or two rows from the same
sequence.

The published algorithm also enforces a **partial order** on blocks: block *A*
precedes *B* if some sequence has *A*'s row before *B*'s, and this relation must
be acyclic. Because `multiz` invocations at a node are independent, cycles do
arise; published TBA restores the order by *deconstructing* offending blocks —
splitting a block into its left-subtree rows and its right-subtree rows. In the
eight-mammal CFTR benchmark this affected 7,134 of 3,548,864 columns, ~0.2%.

**That pass is not active in this distribution.** The only code referencing it,
`iterate_porder()`, sits behind `#ifdef DEAD_CODE` (`tba.c:72`) and calls an
undefined `MD` macro (`tba.c:79-80`) naming a program the Makefile does not
build. Nothing else in the tree mentions partial-order restoration. Treat the
acyclicity property as a description of the published method rather than a
guarantee about this binary's output.

### How the alignment is scored

`multiz` scores a multiple alignment as the sum of all implied pairwise
alignments, with substitutions scored by the HOXD70 matrix
(`mz_scores.c:10-13`):

|   | A | C | G | T |
|---|---:|---:|---:|---:|
| **A** | 91 | −114 | −31 | −123 |
| **C** | −114 | 100 | −125 | −31 |
| **G** | −31 | −125 | 100 | −114 |
| **T** | −123 | −31 | −114 | 91 |

A gap of length *k* costs 400 + 30*k* (`mz_scores.c:23-24`), extended to
multiple alignments with quasi-natural gap costs. Overhangs — gaps at the
extreme ends — are exempt from the **gap-open** charge only (`mz_yama.c:16`,
and the `// not for end-gaps` guards at `mz_yama.c:123,211`); the per-column
extension penalty still applies. The supplement's blanket "not penalized" is
looser than what the code does.

A second matrix calibrated for 85% neutral identity, with gap costs 600 + 50*k*,
is compiled in (`mz_scores.c:26-27`) — but **`init_scores85()` is never
called**. Every program in the suite calls `init_scores70()`, so this matrix is
unreachable here; upstream it belonged to HUMOR's mouse–rat alignments.

**A different scoring matrix cannot be selected per species pair.** If you need
pair-specific scoring, apply it in `lastz` via the specs file, since that is
what determines which regions align at all.

The dynamic programming runs in a band around the guiding alignment's path,
widened by the radius `R`.

## The three steps

1. **Pairwise alignments** — [`all_bz`](programs.md#all_bz) runs `lastz` for
   every species pair, converting and post-processing the results into
   `sp1.sp2.sing.maf` files. It only harvests *names* from the tree argument
   (`all_bz.c:360-378`); the topology is irrelevant to which pairs it runs.
2. **Multiple alignment** — [`tba`](programs.md#tba) or
   [`roast`](programs.md#roast) walks the guide tree leaves-to-root, calling
   `multiz` at each internal node.
3. **Projection** — [`maf_project`](programs.md#maf_project) renders the result
   from the viewpoint of one chosen species.

## Walkthrough: a small TBA alignment

Aligning one locus across human, chimp, mouse, rat, and chicken.

### Set up the working directory

One sequence file per species, named exactly as in the tree, no extension:

```bash
mkdir tba-run && cd tba-run
for s in human chimp mouse rat chicken; do
    ln -s /path/to/data/$s.fa $s
done
```

If any species has more than one FASTA entry, its headers must be
[formatted](formats.md#formatted-headers-multiple-contigs-per-species).

### Write the guide tree

Binary, no branch lengths, spaces for commas, no semicolon:

```bash
TREE="(((human chimp) (mouse rat)) chicken)"
```

Keep every internal node binary. A polytomy such as `((human chimp mouse) ...)`
is **not** an error — it is silently read as `((human chimp) mouse)`, so you
would get an alignment built against a tree you never wrote. See
[the guide tree](formats.md#the-guide-tree).

### Step 1 — pairwise alignments

```bash
all_bz + "$TREE" > all_bz.log 2>&1
```

`+` echoes each command as it runs. This produces `.orig.maf` (raw) and
`.sing.maf` (single-coverage) for all ten species pairs.

On a cluster, generate the commands instead of running them:

```bash
all_bz - "$TREE" lastz.specs > commands.txt
```

then distribute `commands.txt` however your scheduler prefers.

### Step 2 — the multiple alignment

```bash
tba + "$TREE" *.*.sing.maf tba.maf > tba.log 2>&1
```

The sequence files must still be in the directory — `tba` passes their names to
`pair2tb`. The multiz installation directory must be on `PATH`, because `tba`
invokes `multiz`, `maf_project`, `pair2tb`, and `get_covered` by bare name.

### Step 3 — project

```bash
maf_project tba.maf human > human.tba.maf
```

Project onto any other species as needed; the results are guaranteed
consistent with each other:

```bash
for s in chimp mouse rat chicken; do
    maf_project tba.maf $s > $s.tba.maf
done
```

### Then

```bash
# Check the projection for overlapping reference rows
# (this detects overlaps only — it cannot see gaps)
maf_checkThread human.tba.maf

# Drop species and fix row order without re-aligning
maf_order human.tba.maf human mouse rat > human.3way.maf

# Pull out a region
mafFind human.tba.maf 100000 200000 > region.maf

# Convert to FASTA for downstream tools
maf2fasta human human.tba.maf fasta2 > human.tba.fa
```

## tba or roast?

|  | `tba` | `roast` |
|---|---|---|
| Reference | Independent — no species privileged | Required (`E=`) |
| Pairwise alignments needed | *N*(*N*−1)/2 | *N*−1 |
| Output | Threaded blockset, projectable onto any species | Ref-blockset for `E=` only |
| Input alignments preserved | Not guaranteed | Designed to be (Hou 2007, §4.2.4) |
| Runtime | Substantially higher | Substantially lower |

Use **`tba`** when you want a genuine reference-independent blockset and intend
to look at it from more than one species' point of view, and when *N* is small
enough that *N*(*N*−1)/2 pairwise alignments is affordable.

Use **`roast`** when you have a reference in mind and want alignments to it.
This is the genome-scale choice — it is what UCSC runs for their multi-way
conservation tracks.

`roast` guarantees **partial orthology**: every non-reference row in a block is
orthologous to the reference row. It guarantees *complete* orthology — that any
two rows are orthologous — only when the guide tree is **skewed** (each internal
node has a leaf child) and the reference sits at the innermost node.

### Which aligner: multiz or multic?

`multiz` requires that the reference have no duplicates in either input.
`multic` does not. Choose `P=multic` when single coverage has not been enforced
on every species — which is exactly the case under `tba E=species`, where
`all_bz F=species` enforced it on the reference only.

## Driving the pipeline by hand

`all_bz` is a convenience wrapper. Running the steps yourself gives control over
per-pair `lastz` parameters, which is what UCSC does. For one pair:

```bash
blastzWrapper human mouse Y=9000 H=0 $LASTZ_ARGS \
  | lav2maf /dev/stdin human mouse \
  | maf_sort /dev/stdin human \
  > human.mouse.orig.maf

single_cov2 human.mouse.orig.maf > human.mouse.sing.maf
```

This reproduces exactly what `all_bz` generates (`all_bz.c:46,48-49`). Produce
one `.sing.maf` per pair, then run `tba` on the result as usual.

To re-do post-processing without re-aligning — the reason `.orig.maf` files are
kept — use `all_bz b=0`.

## Genome-scale alignment: how UCSC does it

The UCSC Genome Browser's multi-way conservation tracks are built with this
software. Their `multiz.2009-01-21_patched` is the same upstream release this
repository forked, and their build notes in the kent tree
(`src/hg/makeDb/doc/*/multiz*way.txt`) are the most detailed published record of
running it at scale. Two things there are worth knowing.

### The N-way pipeline does not use all_bz

For their multi-way alignments, pairwise input comes from the standard UCSC
chain/net pipeline — `lastz`, then `axtChain`, then `chainNet` — exported as
MAF, renamed to the `db.species.sing.maf` convention, and fed to `autoMZ`.
`all_bz` appears nowhere in the kent tree.

`single_cov2` is skipped too, and the build notes give the reason: "Single-cov
guaranteed by use of net MAF's, so it is not necessary to run single_cov2"
(`src/hg/makeDb/doc/danRer2.txt:4553-4554`). Netting already enforces the
property that `single_cov2` exists to impose.

Which net flavour is used is decided per species by **assembly quality and
distance, in that order** — not by distance alone. The hg38 100-way states the
rule as code (`multiz100way.txt:331-338`):

| Net | Condition | Count in the 100-way |
|---|---|---|
| plain `mafNet` | distance > 0.70 | 42 |
| reciprocal-best | distance ≤ 0.70 **and** poor assembly (>20% N, >2000 contigs, or N50 < 1 Mb) | 8 |
| `synNet` | distance ≤ 0.70 **and** good assembly | 49 |

So the most distant species get the *plainest* net, and reciprocal-best is a
remedy for fragmented assemblies rather than for distance. Smaller builds
simplify this — the *C. elegans* 26-way uses ordinary net for everything
(`ce11/multiz26way.txt:312`).

The naming convention is the whole interface. Anything that can produce a
pairwise MAF can feed this suite.

### The cluster driver

The unit of work is not a chromosome. `mafSplit` cuts the pairwise MAFs into
roughly 5 Mbp chunks at gaps of at least 10 kb
(`multiz30way.txt:363-366`), one `autoMZ` job runs per chunk, and a later pass
reassembles per-chromosome files. The hg38 30-way ran 678 such chunks across
358 chromosomes with data.

Reduced from the real `autoMultiz.csh` (`multiz30way.txt:437-464`), which is
csh; the shape is what matters:

```csh
set db = hg38
set c = $1                        # e.g. hg38_chr1.00.maf — a chunk, not a chrom
set tmp = /dev/shm/$db/multiz.$c
mkdir -p $tmp && cp tree.nh species.list $tmp && cd $tmp

# One .sing.maf per non-reference species, including empty placeholders
foreach s (`sed -e "s/$db //" species.list`)
    set in = $pairs/$s/$c
    set out = $db.$s.sing.maf
    if (-e $in.gz) then
        zcat $in.gz > $out                              # the usual path
        if (! -s $out) echo "##maf version=1 scoring=autoMZ" > $out
    else if (-e $in) then
        ln -s $in $out
    else
        echo "##maf version=1 scoring=autoMZ" > $out
    endif
end

set path = ($run/penn $path); rehash
autoMZ + T=$tmp E=$db "`cat tree.nh`" $db.*.sing.maf $c
```

Four details that generalise:

- **Every species needs a file, even an empty one.** A species with no
  alignments in this chunk still gets a `.sing.maf` containing only the
  `##maf version=1` header. A missing file is an error; an empty one is fine.
  Note the script checks this twice — once for a missing input, once for an
  input that decompresses to nothing.
- **Work on local scratch.** `T=$tmp` points `roast`/`autoMZ` at `/dev/shm`.
  The job copies results out and deletes the directory.
- **`PATH` must contain the multiz binaries**, because the driver shells out to
  `multiz` and `maf_project` by bare name — hence the `set path` line.
- **Budget real time.** The hg38 30-way took 678 cluster jobs, averaging 95
  minutes each, about 1,069 CPU hours, 22 hours wall-clock
  (`multiz30way.txt:484-489`).

### A TBA run at scale

`src/hg/makeDb/doc/hg38/tba10way.txt` records a genuine `tba` run over ten
species, anchored on human chrX. Unlike the N-way pipeline above, this build
*does* use the suite's own pairwise tools: a hand-rolled `runOne` script
(`tba10way.txt:861-889`) runs

```bash
lastzWrapper $target $query $lastzArgs \
  | lav2maf /dev/stdin $target $query \
    | maf_sort /dev/stdin $target > $target.$query.orig.maf
single_cov2 $target.$query.orig.maf > $target.$query.sing.maf
```

for each of **45 all-by-all pairs** (`tba10way.txt:811-812`) — which is exactly
what `all_bz` would have generated. Then:

```bash
tba "(((((((hg38 panTro6) rheMac10) mm10) ((canFam4 neoSch1) pteAle1)) loxAfr3) monDom5) ornAna2)" \
    *.sing.maf chrX.tba10way.maf
# real    154m33.929s  →  1.6 GB of MAF

maf_project chrX.tba10way.maf hg38 > hg38.chrX.tba10way.maf
```

They then project the same blockset onto five other species to obtain mutually
consistent alignments — the property from [Concepts](#concepts), exercised in
production. Note that they copy the FASTA into species-named files in the work
directory (`cp -p fasta/${S}.fa /dev/shm/chrX.10way/${S}`, `tba10way.txt:899`)
before running `tba`, confirming that `tba` needs the sequences and not only
the MAFs.

The binary they invoke is named `lastzWrapper` rather than `blastzWrapper`,
which suggests they made the same lastz substitution this fork does. That is an
inference from the name alone: `lastzWrapper` appears exactly once in the kent
tree and `blastzWrapper` not at all, and no source or install line for it is
published.

## Troubleshooting

**`roast` dies with a `multiz` usage dump** — you probably passed `C=` without
`P=multic`. `multiz` does not accept `C=`, and rather than naming the flag it
falls through to an argument-count error. Note `C=` has no effect even with
`P=multic` in this build; see [`multic`](programs.md#multic).

**A sub-command reports a missing file, but the file you expect is right
there** — check the name. `tba` builds `sp1.sp2.sing.maf` from the tree, in the
tree's order. Also check that the *sequence* files exist, named by species with
no extension.

**`sh: multiz: command not found` from inside `tba`** — the installation
directory must be on `PATH`. Having the binaries in the current directory is
not enough.

**`single_cov2` produced no deleted-regions file** — `F=` must come *after*
`R=`. Given in the other order it is silently ignored. See
[`single_cov2`](programs.md#single_cov2).

**Specs-file options appear to have no effect** — the parser is particular about
tabs versus spaces and reports nothing when it misreads the file. Check the
whitespace on the indented option lines.

**`maf_checkThread` reports no errors on an alignment you suspect is broken** —
it only detects *overlaps*, never gaps, and it checks the first row only. A
projection with holes in the reference passes silently. See
[`maf_checkThread`](programs.md#maf_checkthread).

**Blocks you expected are missing** — `maf_order` drops single-row blocks
unless you pass `all`. `multiz` and `multic` accept `all` too, but it does
nothing there: they always emit single-row blocks. See
[`multiz`](programs.md#multiz).

**`toast: command not found`** — you selected `all_bz A=0` or `A=2`. Those
paths need the separate TOAST distribution; see
[External programs](programs.md#external-programs-not-included).

**`tba`: `no alignment found for X and Y` after setting `X=1` or `X=2`** — those
options do not run `toast`; they only change which *suffix* `tba` looks for
(`.toast.maf`, `.toast2.maf`). If you have not produced those files, `tba`
simply cannot find its inputs.

## References

- Blanchette M, Kent WJ, Riemer C, Elnitski L, Smit AFA, Roskin KM, Baertsch R,
  Rosenbloom K, Clawson H, Green ED, Haussler D, Miller W. "Aligning multiple
  genomic sequences with the threaded blockset aligner." *Genome Research*
  14:708–715, 2004. doi:10.1101/gr.1933104
- Margulies EH, Hou M, Miller W. "A Practical Guide to Using TBA."
  April 28, 2005.
- Hou M. "Algorithms for aligning and clustering genomic sequences that contain
  duplications." PhD thesis, Penn State, August 2007. — source for the `roast`
  algorithm and the TOAST/BOAST context.
- UCSC kent source tree, `src/hg/makeDb/doc/*/multiz*way.txt` and
  `src/hg/makeDb/doc/hg38/tba10way.txt`.
