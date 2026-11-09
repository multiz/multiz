#!/usr/bin/env python3
"""Five seeded synthetic species for the roast and tba tests: one 10 kb ancestor and its
descendants on the tree ((sa sb) (sc (sd se))), with substitutions, small indels, two 2 kb
inversions (long enough that lastz keeps them as - strand blocks), a short inversion, a moved
segment, a duplication and a novel insert. Usage: make_species.py OUTDIR (one FASTA per species)."""
import random
import sys

rnd = random.Random(20261107)
comp = str.maketrans("ACGT", "TGCA")


def mutate(s, sub, indel):
    out = []
    for ch in s:
        r = rnd.random()
        if r < indel / 2:
            continue                                            # deletion
        if r < indel:
            out.append(rnd.choice("ACGT") * rnd.randint(1, 6))  # insertion
        out.append(rnd.choice("ACGT") if rnd.random() < sub else ch)
    return "".join(out)


def invert(s, i, j):
    return s[:i] + s[i:j].translate(comp)[::-1] + s[j:]


anc = "".join(rnd.choice("ACGT") for _ in range(10_000))
ab = mutate(anc, 0.03, 0.004)
cde = mutate(anc, 0.05, 0.006)
de = mutate(cde, 0.03, 0.004)
seqs = {"sa": mutate(ab, 0.02, 0.003), "sb": mutate(ab, 0.02, 0.003), "sc": mutate(cde, 0.04, 0.004),
        "sd": mutate(de, 0.02, 0.003), "se": mutate(de, 0.02, 0.003)}
seqs["sd"] = invert(invert(seqs["sd"], 3_000, 5_000), 7_500, 7_625)
seqs["sb"] = invert(seqs["sb"], 1_500, 3_500)
s = seqs["sc"]
seqs["sc"] = s[:1_000] + s[6_000:6_750] + s[1_000:6_000] + s[6_750:]     # a moved segment
s = seqs["se"]
seqs["se"] = s[:3_000] + s[2_000:2_400] + s[3_000:]                      # a duplication
s = seqs["sa"]
seqs["sa"] = s[:8_500] + "".join(rnd.choice("ACGT") for _ in range(200)) + s[8_500:]  # a novel insert
for name, s in seqs.items():
    with open(f"{sys.argv[1]}/{name}", "w") as fh:
        fh.write(f">{name}\n")
        for k in range(0, len(s), 70):
            fh.write(s[k:k + 70] + "\n")
