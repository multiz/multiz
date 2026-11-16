"""Write random yama cases (see main.c) to stdout: gen.py SEED NCASES [MAXDIM] [MODE]

MODE random (default): columns mix ACGT, lower case, N and dash runs (including all-dash
columns); about a quarter of cases use the full band, the rest a random monotonic band of
width at least MIN(N,10), which exercises yama's LB/RB edge guards.

Targeted modes (tests/yama_fuzz/check.sh), aimed at what random cases reach rarely:
  narrow  every band width the vector loops can see at the minimum: N = 1..9 with the full
          band (width N), and constant widths 10..17 (every tail length mod 4), with LB steps
          of 0-3 so cells next to the band edge read MININT sentinels
  ties    columns from a single base, or all N, with M != N: many equal-scoring paths, so
          exact C/D/I ties, including D/I ties (the tie_di mutant), decide the traceback
  edges   M or N in 1..3, all-dash columns at the first and last positions, and B starting
          with dash-heavy columns: rows 0 and M and the first-column gap-open waiver
"""
import random, sys

def column(rng, h, dashy):
    return "".join("-" if rng.random() < dashy else rng.choice("ACGTACGTacgtN") for _ in range(h))

def block(rng, h, n):
    dashy = rng.choice([0.0, 0.1, 0.3, 0.6])
    cols, run = [], 0
    for _ in range(n):
        if run == 0 and rng.random() < 0.1:
            run = rng.randint(1, 8)
        if run:
            run -= 1
            c = list(column(rng, h, dashy))
            for i in rng.sample(range(h), rng.randint(1, h)):
                c[i] = "-"
            cols.append("".join(c))
        else:
            cols.append(column(rng, h, dashy))
    return cols

def band(rng, M, N):
    mw = min(N, 10)
    if rng.random() < 0.25:
        return [0] * (M + 1), [N] * (M + 1)
    w = rng.randint(mw, max(mw, min(N, rng.choice([10, 20, 40, 70, 150]))))
    LB, RB = [0], [min(N, rng.randint(mw, max(mw, w)))]
    for r in range(1, M + 1):
        lb = min(LB[-1] + rng.choice([0, 0, 1, 1, 1, 2, 3]), N - mw, RB[-1])
        rb = min(N, max(RB[-1], lb + mw, lb + w + rng.randint(-3, 3)))
        LB.append(lb); RB.append(rb)
    RB[M] = N
    return LB, RB

def narrow_band(rng, M, N):
    if N <= 9:
        return [0] * (M + 1), [N] * (M + 1)
    w = rng.randint(10, min(17, N))
    LB, RB = [0], [w]
    for r in range(1, M + 1):
        lb = min(LB[-1] + rng.choice([0, 0, 1, 1, 2, 3]), N - w, RB[-1])
        LB.append(lb); RB.append(min(N, max(RB[-1], lb + w)))
    RB[M] = N
    return LB, RB

def ties_block(rng, h, n):
    kind = rng.choice(["A", "N", "Aa", "A-"])
    if kind == "A-":   # one base, with whole-dash columns sprinkled in
        return ["-" * h if rng.random() < 0.2 else "A" * h for _ in range(n)]
    c = {"A": "A", "N": "N", "Aa": "a"}[kind]
    return [c * h for _ in range(n)]

def edges_block(rng, h, n):
    cols = block(rng, h, n)
    if n >= 1 and rng.random() < 0.5: cols[0] = "-" * h
    if n >= 2 and rng.random() < 0.5: cols[-1] = "-" * h
    return cols

def targeted(rng, mode, out):
    K, L = rng.randint(1, 4), rng.randint(1, 4)
    if mode == "narrow":
        M, N = rng.randint(1, 60), rng.choice([rng.randint(1, 9), rng.randint(10, 80)])
        LB, RB = narrow_band(rng, M, N)
        A, B = block(rng, K, M), block(rng, L, N)
    elif mode == "ties":
        M = rng.randint(1, 40); N = rng.randint(1, 40)
        if N == M: N += rng.choice([1, 2, 5])
        LB, RB = band(rng, M, N) if rng.random() < 0.5 else ([0] * (M + 1), [N] * (M + 1))
        A, B = ties_block(rng, K, M), ties_block(rng, L, N)
    else:  # edges
        M, N = rng.choice([(rng.randint(1, 3), rng.randint(1, 40)), (rng.randint(1, 40), rng.randint(1, 3)),
                           (rng.randint(1, 3), rng.randint(1, 3))])
        LB, RB = [0] * (M + 1), [N] * (M + 1)
        A = edges_block(rng, K, M)
        B = edges_block(rng, L, N)
        if N >= 2 and rng.random() < 0.5:   # B opens with dash-heavy columns: the first-column waiver
            B[0] = "-" * (L - 1) + "A" if L > 1 else "-"
    out.append(f"{K} {M}"); out += A
    out.append(f"{L} {N}"); out += B
    out.append(" ".join(map(str, LB))); out.append(" ".join(map(str, RB)))

def main():
    seed, ncases = int(sys.argv[1]), int(sys.argv[2])
    maxdim = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    mode = sys.argv[4] if len(sys.argv) > 4 else "random"
    rng = random.Random(seed)
    out = []
    if mode != "random":
        if mode not in ("narrow", "ties", "edges"):
            sys.exit(f"gen.py: unknown mode {mode}")
        for _ in range(ncases):
            targeted(rng, mode, out)
        sys.stdout.write("\n".join(out) + "\n")
        return
    for _ in range(ncases):
        K, L = rng.randint(1, 8), rng.randint(1, 8)
        if rng.random() < 0.1:
            K, L = rng.randint(1, 60), rng.randint(1, 60)
        M, N = rng.randint(1, maxdim), rng.randint(1, maxdim)
        if rng.random() < 0.3:
            M, N = rng.randint(1, 12), rng.randint(1, 12)
        LB, RB = band(rng, M, N)
        out.append(f"{K} {M}"); out += block(rng, K, M)
        out.append(f"{L} {N}"); out += block(rng, L, N)
        out.append(" ".join(map(str, LB))); out.append(" ".join(map(str, RB)))
    sys.stdout.write("\n".join(out) + "\n")

main()
