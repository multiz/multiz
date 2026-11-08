#!/usr/bin/env python3
"""Write a pair of synthetic MAF files for multiz: one with K rows per block, one with L, sharing a reference row.

    python3 make_blocks.py OUT_PREFIX K L [--blocks 10] [--width 2000] [--seed 1]

Writes OUT_PREFIX.1.maf (the reference plus K - 1 species) and OUT_PREFIX.2.maf (the reference plus L - 1
others). Each species row is the reference with 2% substitutions and 1% single-base deletions. Each file also
gets its own insertion columns (a dash in the reference row), at different places in the two files. So the two
files disagree about the columns, and multiz must align them with its dynamic programming (`yama`). The output
is deterministic for a given seed.
"""
from __future__ import annotations

import argparse
import random

BASES = "ACGT"


def species_row(ref: str, rng: random.Random) -> str:
    out = []
    for c in ref:
        r = rng.random()
        if r < 0.01:
            out.append("-")
        elif r < 0.03:
            out.append(rng.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def add_insertions(rows: list[str], rng: random.Random, rate: float = 0.01) -> list[str]:
    """Insert columns where the reference has a dash and the other rows have a base or a dash."""
    cols = list(zip(*rows))
    out = []
    for col in cols:
        out.append(col)
        if rng.random() < rate:
            out.append(("-",) + tuple(rng.choice(BASES + "-") for _ in col[1:]))
    return ["".join(r) for r in zip(*out)]


def write(path: str, names: list[str], blocks: list[tuple[int, list[str]]], sizes: dict[str, int]) -> None:
    pos = {n: 0 for n in names}
    with open(path, "w") as fh:
        fh.write("##maf version=1\n")
        for start, rows in blocks:
            fh.write("a score=0\n")
            for i, (name, text) in enumerate(zip(names, rows)):
                n = len(text) - text.count("-")
                s = start if i == 0 else pos[name]
                fh.write(f"s {name} {s} {n} + {sizes[name]} {text}\n")
                pos[name] = s + n
            fh.write("\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("prefix")
    ap.add_argument("K", type=int)
    ap.add_argument("L", type=int)
    ap.add_argument("--blocks", type=int, default=10)
    ap.add_argument("--width", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    gap = 100
    total = a.blocks * (a.width + gap)
    names1 = ["ref.chr1"] + [f"a{i}.chr1" for i in range(1, a.K)]
    names2 = ["ref.chr1"] + [f"b{i}.chr1" for i in range(1, a.L)]
    sizes = {n: 10 * total for n in names1 + names2}
    sizes["ref.chr1"] = total
    b1, b2 = [], []
    for k in range(a.blocks):
        ref = "".join(rng.choice(BASES) for _ in range(a.width))
        start = k * (a.width + gap)
        r1 = add_insertions([ref] + [species_row(ref, rng) for _ in names1[1:]], rng)
        r2 = add_insertions([ref] + [species_row(ref, rng) for _ in names2[1:]], rng)
        b1.append((start, r1))
        b2.append((start, r2))
    write(f"{a.prefix}.1.maf", names1, b1, sizes)
    write(f"{a.prefix}.2.maf", names2, b2, sizes)


if __name__ == "__main__":
    main()
