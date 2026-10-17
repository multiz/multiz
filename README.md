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
18 executables in the source directory.

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

The Makefile compiles at `-O0` with `-Wall -Wextra`. If you are aligning
anything large, raise the optimisation level — these are
compute-bound programs:

```bash
make CFLAGS="-Wall -Wextra -O2"
```

## Dependencies

At runtime:

- **[lastz](https://github.com/lastz/lastz)** must be on `PATH`. It is invoked
  by `blastzWrapper`, and therefore by `all_bz`.
- **The installation directory must be on `PATH`.** `tba`, `roast`, and
  `all_bz` invoke other programs in the suite by bare name — see the
  [per-driver table](docs/programs.md#external-programs-not-included) — along
  with `cat`, `cp`, `echo`, `grep`, `mv`, `rm`, and `touch`.

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

This is a functionally equivalent copy of the multiz/TBA tarball posted at the
Miller Lab website on **21 January 2009**
(`http://www.bx.psu.edu/miller_lab`). UCSC runs that same release in
production, installed as `multiz.2009-01-21_patched`. Their `_patched` suffix
marks local modifications. The patch is unpublished, but has been reconstructed
and verified against their binary: the substantive change is in `multiz`'s
unused-block handling — see [`multiz`](docs/programs.md#multiz). No later
upstream version appears anywhere in their build notes.

Changes made here:

- **`lastz` replaces `blastz`.** `blastzWrapper` now executes `lastz`
  (`blastzWrapper.c:14`). Blastz has been deprecated for well over a decade and
  lastz is a drop-in replacement. The program keeps its original name so that
  existing scripts and the rest of the suite continue to work. Lastz accepts the
  blastz option letters this suite passes (`Y=`, `H=`, `Q=`), which is what
  makes the substitution work. UCSC appears to have done the same — their TBA
  build notes invoke a `lastzWrapper` — though they publish no source for it.
- Syntax corrections to silence compiler warnings; removal of obsolete `rcsid`
  tags. No behavioural change.
- Addition of the MIT license and this documentation.

The motivation for the repository is to provide an equivalent of that tarball
that modern package managers can install.

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
