# Input and output formats

Everything in this suite is glued together by filename conventions and by two
file formats: FASTA for sequence, MAF for alignments. The programs do very
little validation, so a malformed header or a misnamed file usually produces a
confusing error or silently wrong output rather than a clear complaint. This
page is the reference for getting those details right.

- [Species names and file naming](#species-names-and-file-naming)
- [Sequence files (FASTA)](#sequence-files-fasta)
- [The guide tree](#the-guide-tree)
- [The lastz specs file](#the-lastz-specs-file)
- [Alignment files (MAF)](#alignment-files-maf)
- [Pairwise MAF naming conventions](#pairwise-maf-naming-conventions)

## Species names and file naming

A *species name* is the single identifier that ties together every part of a
run. For a species named `human`:

| Thing | Must be named |
|---|---|
| Sequence file | `human` (no extension) |
| Name in the guide tree | `human` |
| Name in the lastz specs file | `human` |
| `string1` in formatted FASTA headers | `human` |
| Pairwise alignment against `mouse` | `human.mouse.sing.maf` |

`tba` and `roast` build these filenames by string concatenation. A missing
*pairwise MAF* is caught by `tba` itself — `no alignment found for X and Y`
(`tba.c:125-127`). A missing or misnamed *sequence file* is not: the name is
passed to `pair2tb` and surfaces as an error from that child process instead.
(`roast` never calls `pair2tb`; it uses only `multiz`/`multic` and
`maf_project`.)

Two practical consequences:

- **Sequence files need no extension.** If your data lives in `human.fa`,
  create a link: `ln -s human.fa human`. This is exactly what UCSC does in
  their production runs.
- **Keep names short, alphabetic-initial, and free of `(`, `)`, and
  whitespace.** `tba` and `roast` parse the tree with `speciesTree.c`, which
  requires a name to *begin* with a letter (`speciesTree.c:54`), allows letters,
  digits, `_` and `.` thereafter (`speciesTree.c:59`), and treats any other
  character as fatal (`speciesTree.c:69-70`). `all_bz` uses a looser tokeniser
  of its own (`NON_NAME " ()"`, `all_bz.c:44`) — so a name starting with a digit
  passes `all_bz` and then fails in `tba`. Avoid `.` as well: it is the
  separator in generated MAF filenames.

## Sequence files (FASTA)

One file per species, containing one or more FASTA entries.

### Simple case: one sequence per species

If a species' file contains exactly one FASTA entry, the header content does
not matter and need not follow any convention. This is the common case for
aligning a single locus across species.

```
>human
ACGTACGTACGT...
```

### Formatted headers: multiple contigs per species

If a species contributes more than one FASTA entry — unassembled contigs, a
syntenic break, several chromosomes — every entry must carry a header that
uniquely identifies it:

Two forms are accepted. The parser tries the **range** form first
(`multi_util.c:315`), then falls back to the **start-only** form
(`multi_util.c:317`):

```
>string1:string2:start-end:char:int2      <- range form, tried first
>string1:string2:int1:char:int2           <- start-only form
```

| Field | Meaning |
|---|---|
| `string1` | Source organism identifier. **Must equal the sequence file's name.** |
| `string2` | Contig/chromosome name. Must be unique within the species' file. |
| `int1` / `start` | Start offset of this sequence in the full contig, 1-based and inclusive. Use `1` if the entry is the whole contig. |
| `end` | Range form only. End offset, 1-based and inclusive. In the start-only form this is computed as `start + length - 1`. |
| `char` | Strand, `+` or `-`. Almost always `+`. |
| `int2` | Size of the entire contig named by `string2` — not the length of this entry. |

Note that `parse_fasta` (`multi_util.c:89-91`) accepts **only** the range form,
so prefer it if you hit an unexplained header error.

`string1.string2` becomes the `src` field of the MAF components, so a header of
`>human:chr7:1:+:158821424` yields components named `human.chr7`.

The `int1`/`int2` fields exist so that a sub-sequence of a chromosome carries
genome coordinates through the alignment. If you never work with
sub-sequences, set `int1` to `1` and `int2` to the sequence's own length.

`get_standard_headers` reports the `start-end:strand:srcSize` portion for each
entry of a plain FASTA file, assuming each entry starts at position 1. It emits
`1-<len>:+:<len>` — a hyphen between start and end, and `end` equal to
`srcSize`. It is a report, not a ready-made header; you still compose the
`species:contig:` prefix yourself. See
[`get_standard_headers`](programs.md#get_standard_headers).

### ENCODE MSA headers

The pipe-delimited header format adopted by the ENCODE MSA group is also
accepted:

```
>${COMMON_NAME}|${ENCODE_REGION}|${FREEZE_DATE}|${NCBI_TAXON_ID}|${ASSEMBLY_PROVIDER}|${ASSEMBLY_DATE}|${ASSEMBLY_ID}|${CHROMOSOME}|${CHROMOSOME_START}|${CHROMOSOME_END}|${CHROM_LENGTH}|${STRAND}|${ACCESSION}.${VERSION}|${NUM_BASES}|${NUM_N}|${THIS_CONTIG_NUM}|${TOTAL_NUM_CONTIGS}|${OTHER_COMMENTS}
```

Empty fields are represented by `.`.

## The guide tree

A binary tree in a **modified Newick format**: branch lengths removed, commas
replaced by spaces, trailing semicolon removed, and the whole thing passed as a
single shell argument.

```
"(((((((human chimp) gorilla) baboon) (rat mouse)) (cow pig)) chicken) fugu)"
```

Requirements:

- **Make it strictly binary.** Every internal node should have exactly two
  children. A polytomy is **not rejected** — it is silently re-interpreted.
  `speciesTree.c:71` merges the top two completed subtrees after every token,
  so `(a b c)` is quietly read as `((a b) c)`, a left comb. You get an
  alignment, computed against a tree you did not specify.
- **Leaf names are species names**, matching sequence filenames exactly.
- **No branch lengths, no commas, no semicolon.**
- Quote it, because it contains spaces and parentheses.

To convert a standard Newick tree with branch lengths, UCSC uses:

```bash
sed 's/[a-z][a-z]*_//g; s/:[0-9\.][0-9\.]*//g; s/;//; /^ *$/d' input.nh \
  | xargs echo | sed 's/ //g; s/,/ /g' > tree.nh
```

and derives the species list from the result with:

```bash
sed 's/[()]//g; s/,/ /g' tree.nh > species.list
```

The tree's *shape* matters for `roast`: only when the tree is skewed (a
caterpillar, each internal node having a leaf child) and the reference sits at
the innermost node are all pairs of rows in a block guaranteed orthologous.
See [`docs/pipelines.md`](pipelines.md#tba-or-roast).

## The lastz specs file

Optional. Supplies per-species-pair command-line options to `lastz`. Referred
to throughout the original documentation as "the blastz specs file".

```
# This is a sample specs file
#define PLACENTAL human chimp gorilla baboon rat mouse cow pig
#define NON_PLACENTAL chicken fugu
PLACENTAL : PLACENTAL
	B=2 C=0
NON_PLACENTAL : *
	Q=HoxD55.q
```

- `#define NAME members...` declares a named group of species.
- `GROUP : GROUP` introduces a rule; `*` matches any species.
- The **indented** line following a rule holds the options passed to `lastz`
  for pairs matching that rule.

There is a real whitespace trap, but it is on the **`#define` line**, not the
indented one. `all_bz.c:106-108` searches the whole line for a space *before*
looking for a tab, so `#define PLACENTAL<TAB>human chimp` splits at the space
between `human` and `chimp` — yielding a macro named `PLACENTAL\thuman` with
the single member `chimp`, silently wrong. **Separate the macro name from its
members with a space, not a tab.**

The indented option line is checked and does complain: a line that does not
begin with whitespace is a fatal `missing space at start of ...`
(`all_bz.c:143-147`), and there either a tab or spaces will do.

Options given here are appended to the ones `all_bz` always supplies
(`Y=9000 H=0`, see `all_bz.c:46`).

## Alignment files (MAF)

Multiple Alignment Format, as specified by UCSC:
<https://genome.ucsc.edu/FAQ/FAQformat.html#format5>

A file is a `##maf` header line followed by blocks. Each block is an `a` line
carrying a score, then one `s` line per row:

```
##maf version=12 scoring=tba.v12
a score=23262.0
s human.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA...
s mouse.chr6    53215344 38 + 151104725 -AAAGGGAATGTTAAGCAAACGA...
s rat.chr4      81344243 40 + 187371129 -AAAGGGGATGCTAAGCCAATGA...
```

Fields on an `s` line: `src`, `start`, `size`, `strand`, `srcSize`, and the
gapped sequence text.

The `version=` number above is not a typo: `tba` substitutes its own program
version into that field (`tba.c:409`), so it writes `version=12` where the MAF
spec expects `version=1`. `roast` writes `version=1` correctly
(`auto_mz.c:266`). Parsers in this suite do not care; stricter third-party
parsers might.

Two conventions that regularly trip people up:

- **MAF coordinates are 0-based**, while the `int1` field of a formatted FASTA
  header is 1-based and inclusive. The conversion happens when the alignment is
  built.
- **`start` is always relative to the given strand.** For a `-` strand row,
  `start` counts from the end of the source sequence.

A block whose rows all come from distinct species, with the reference first, is
what these programs call a *block* of a *blockset*. A block may have a single
row — such blocks represent sequence present in the blockset but aligned to
nothing. `multiz` and `multic` always emit them (their `all` flag is inert in
this build); `maf_order` drops them unless you pass `all`.

## Pairwise MAF naming conventions

`all_bz` writes, and `tba`/`roast` read, files named `species1.species2.SUFFIX`
where the suffix records how far post-processing has gone:

| Suffix | Produced by | Selected by |
|---|---|---|
| `.orig.maf` | `lastz` → `lav2maf` → `maf_sort` | — (always written) |
| `.sing.maf` | `single_cov2` on the `.orig.maf` | `X=0` (default) |
| `.clean.maf` | `blastz_clean` (**not included**) | — |
| `.toast.maf` | `toast` (**not included**) | `X=1` |
| `.toast2.maf` | `chain` (**not included**) | `X=2` |

The `.orig.maf` intermediates are kept deliberately: re-running post-processing
with different settings is much cheaper than re-running the aligner, and
`all_bz b=0` does exactly that.

`toast`, `chain`, and `blastz_clean` are **not part of this package**. The
`X=1` and `X=2` paths, and `all_bz A=0` / `A=2`, require the separate TOAST
distribution. See [`docs/programs.md`](programs.md#external-programs-not-included).

The species order in the filename matters: the first named species is the one
whose rows come first in each block. `tba` derives the expected name from the
guide tree, so let `all_bz` generate the files rather than naming them by hand.
