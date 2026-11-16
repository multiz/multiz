"""Machine proof of i_avx2 (mz_yama.c) against the scalar pass 2, with z3.

Usage: python3 scan_z3.py [--self-test] [--equiv] [LEMMA ...]
  LEMMA is any of L0 L1 L2 L3 L4 (default: all).  Runs in under a second.

The proof is in two layers.

Intrinsics, on bit-vectors at their own granularity (bytes for blendv_epi8,
32-bit elements for blend_epi32, 64-bit lanes for permute4x64/cmpgt, the low
signed 32 bits for mul_epi32), for ALL inputs:
  L0a max_epi64 (blendv_epi8 on a cmpgt_epi64 mask) is the lane-wise max.
  L0b blend_epi32 with imm 0x03 / 0x0F moves whole lanes 0 / 0-1.
  L0c permute4x64 with imm 0x90 / 0x40 / 0xFF shifts by one lane / two lanes
      / broadcasts lane 3.
  L1  BEST3_V equals BEST3: value, 2-bit flag, and all-0/all-1 lane masks.
  L2  mul_epi32 equals the scalar 64-bit product when the count (a long long
      copy of an int count) is in [0, INT_MAX] and the row coefficient
      passes FITS32 (L2a-c); L2n is a witness that it differs without FITS32.

The scan, as a lane-level model (scan_lanes) that follows i_avx2's loop body
line for line, with the intrinsics replaced by their L0/L1 meanings:
  L3  In exact integers, every lane's I, the I(col-1) handed to BEST3_V, the
      stored value and the new carry equal the serial recurrence
      I(c) = max(max(x, y), I(c-1) - Z) - E, and so does each scalar step.
  L4  Every add and sub of the block and of 4 scalar steps, evaluated
      exactly, lies in [LLONG_MIN, LLONG_MAX].  By induction over the
      operations in program order, the 64-bit code then sees exact operands,
      returns exact results, and its signed compares agree with the integer
      ones, so L3 holds for the 64-bit vector and scalar code.
  L3n With values allowed below the NEG fill, L3 fails, so the bounds matter.

Bounds for L3/L4, stated rather than derived, because no range check
enforces them: x, y and the carried I in [VLO, VHI] = [-2^62 + 2^47,
2^62 - 2^47]; Z (Izc, a score_t) in [-2^44, 2^44]; E (eK, the int product
n*K*gap_extend) in [INT_MIN, INT_MAX].  MININT = -2^61, so a value would have
to sink 2^61 - 2^47 below MININT to leave them.  YAMA_SHADOW records the
ranges actually reached (on piece 03: values down to MININT - 525,200,
|Izc| <= 471,200, |eK| <= 37,800).

Not machine-checked: that scan_lanes mirrors the C source.  An encoding of the
whole loop body at the intrinsic level, compared with scan_lanes for all
inputs, did not finish within 25 minutes, so the correspondence rests on
L0-L2 plus reading.  Shadow mode and the mutation tests cover it
empirically.

--self-test: three mis-encodings of the block (wrong permute, NEG fill
replaced by ZERO, carry dropped) must each be refuted.

--equiv: the mutant that blends only the low 32 bits of NEG in step 1
(imm 0x01 for 0x03) is equivalent when Z >= 0 and E >= 0 (Izc is a sum of
counts times K*GAP, GAP is 0 or gap_open >= 0, and eK is n*K*gap_extend),
and not equivalent when they may be negative.
"""
import sys
import time
from z3 import (is_bv, Int, IntVal, BitVec, BitVecVal, Concat, Extract, If, SignExt, Solver, And,
                Or, Not, simplify, sat, unsat)

LLONG_MIN = -(1 << 63)
INT_MIN, INT_MAX = -(1 << 31), (1 << 31) - 1
VLO, VHI = -(1 << 62) + (1 << 47), (1 << 62) - (1 << 47)
ZLO, ZHI = -(1 << 44), 1 << 44


def lane(v, i):
    return Extract(64 * i + 63, 64 * i, v)


def pack(lanes):  # lanes[0] is the low lane
    return Concat(*reversed(lanes))


def set1(x):
    return pack([x] * 4)


def lanewise(f, *vs):
    return pack([f(*[lane(v, i) for v in vs]) for i in range(4)])


def add64(a, b):
    return lanewise(lambda x, y: x + y, a, b)


def sub64(a, b):
    return lanewise(lambda x, y: x - y, a, b)


def cmpgt64(a, b):
    return lanewise(lambda x, y: If(x > y, BitVecVal(-1, 64), BitVecVal(0, 64)), a, b)


def bor(a, b):
    return a | b


def band(a, b):
    return a & b


def bandnot(a, b):  # _mm256_andnot_si256(a, b) = ~a & b
    return ~a & b


def blendv8(a, b, m):  # byte j from b when bit 7 of m's byte j is set
    out = []
    for j in range(32):
        out.append(If(Extract(8 * j + 7, 8 * j + 7, m) == 1,
                      Extract(8 * j + 7, 8 * j, b), Extract(8 * j + 7, 8 * j, a)))
    return Concat(*reversed(out))


def blend32(a, b, imm):  # 32-bit element j from b when bit j of imm is set
    out = []
    for j in range(8):
        src = b if (imm >> j) & 1 else a
        out.append(Extract(32 * j + 31, 32 * j, src))
    return Concat(*reversed(out))


def permute4x64(a, imm):
    return pack([lane(a, (imm >> (2 * i)) & 3) for i in range(4)])


def mul_epi32(a, b):  # signed low 32 bits of each lane, exact 64-bit product
    return lanewise(lambda x, y: SignExt(32, Extract(31, 0, x)) * SignExt(32, Extract(31, 0, y)), a, b)


def max_epi64(a, b):
    return blendv8(a, b, cmpgt64(b, a))


def best3_v(x, y, z):
    bx = bor(cmpgt64(y, x), cmpgt64(z, x))
    md = band(bx, cmpgt64(y, z))
    mi = bandnot(md, bx)
    V = blendv8(blendv8(x, y, md), z, mi)
    return V, mi, md


def best3(x, y, z):  # the scalar macro: value and flag (C=0, I=1, D=2)
    bx = And(x >= y, x >= z)
    by = And(Not(bx), y > z)
    V = If(bx, x, If(by, y, z))
    F = If(bx, BitVecVal(0, 2), If(by, BitVecVal(2, 2), BitVecVal(1, 2)))
    return V, F


def prove(name, claim, assumptions=(), expect=unsat):
    t0 = time.time()
    s = Solver()
    s.add(*assumptions)
    s.add(Not(claim))
    r = s.check()
    ok = r == expect
    tag = ("PROVED " if expect == unsat else "REFUTED") if ok else "FAILED "
    print(f"{tag} {name}  [{r}, {time.time() - t0:.1f} s]", flush=True)
    if r == sat:
        m = s.model()
        print("  model:", sorted(((d.name(), m[d].as_signed_long() if hasattr(m[d], "as_signed_long") else m[d])
                                 for d in m.decls()), key=lambda kv: kv[0])[:12], flush=True)
    return ok


def lemma_best3():
    x, y, z = (BitVec(n, 256) for n in "xyz")
    V, mi, md = best3_v(x, y, z)
    cl = []
    for i in range(4):
        sv, sf = best3(lane(x, i), lane(y, i), lane(z, i))
        fbit_i = Extract(63, 63, lane(mi, i))
        fbit_d = Extract(63, 63, lane(md, i))
        cl += [lane(V, i) == sv, Concat(fbit_d, fbit_i) == sf,
               Or(lane(mi, i) == 0, lane(mi, i) == BitVecVal(-1, 64)),
               Or(lane(md, i) == 0, lane(md, i) == BitVecVal(-1, 64))]
    return prove("L1 BEST3_V == BEST3 (value, flags, full-lane masks), all inputs", And(*cl))


def lemma_xy():
    g, c = BitVec("count", 64), BitVec("coef", 64)
    ge, ce = SignExt(32, Extract(31, 0, g)), SignExt(32, Extract(31, 0, c))
    ok = prove("L2a count in [0,INT_MAX]: SignExt(low32(count)) == count", ge == g, [g >= 0, g <= INT_MAX])
    ok &= prove("L2b coef FITS32: SignExt(low32(coef)) == coef", ce == c, [c >= INT_MIN, c <= INT_MAX])
    # L2a and L2b make mul_epi32's operands the scalar operands, so the products
    # are the same term; z3 checks that step by congruence.
    ok &= prove("L2c given L2a, L2b: mul_epi32(count, coef) == count*coef", ge * ce == g * c, [ge == g, ce == c])
    w = simplify(SignExt(32, Extract(31, 0, BitVecVal(1, 64))) * SignExt(32, Extract(31, 0, BitVecVal(1 << 31, 64))))
    good = w.as_signed_long() == -(1 << 31)
    print(f"{'REFUTED' if good else 'FAILED '} L2n without FITS32: count 1, coef 2^31 gives "
          f"{w.as_signed_long()} in the vector code, 2147483648 in the scalar code", flush=True)
    return ok and good


def lemma_intrinsics():
    a, b = BitVec("a", 256), BitVec("b", 256)
    ok = prove("L0a max_epi64 (blendv_epi8 on a cmpgt mask) == lane-wise signed max, all inputs",
               max_epi64(a, b) == lanewise(lambda p, q: If(q > p, q, p), a, b))
    for imm, n in ((0x03, 1), (0x0F, 2)):
        ok &= prove(f"L0b blend_epi32 imm {imm:#04x} == whole lanes 0..{n - 1} from the second operand",
                    blend32(a, b, imm) == pack([lane(b, i) if i < n else lane(a, i) for i in range(4)]))
    for imm, src in ((0x90, (0, 0, 1, 2)), (0x40, (0, 0, 0, 1)), (0xFF, (3, 3, 3, 3))):
        ok &= prove(f"L0c permute4x64 imm {imm:#04x} == lanes {src}",
                    permute4x64(a, imm) == pack([lane(a, j) for j in src]))
    return ok


def scan_lanes(x, y, Z, E, carry, mutant=None, const=lambda v: BitVecVal(v, 64)):
    """i_avx2's loop body after x and y, on 64-bit lanes (index 0 = low lane),
    using L0a-L0c to replace max_epi64, blend_epi32 and permute4x64 by their
    lane-level meaning, and L1 for BEST3_V.  Returns (a, Ip, V, carry, ops)."""
    NEG = const(LLONG_MIN // 2)
    Z0 = const(0)
    ops = []

    def mx(p, q):
        return If(q > p, q, p)

    def add(p, q):
        ops.append(("add", p, q))
        return p + q

    def sub(p, q):
        ops.append(("sub", p, q))
        return p - q

    def shift(v, k, fill):  # permute (0x90: k=1, 0x40: k=2) then blend fill into the low k lanes
        return [fill(i, v) if i < k else v[i - k] for i in range(4)]

    neg = (lambda i, v: Z0) if mutant == "negzero" else (lambda i, v: NEG)
    if mutant == "blend01":  # only the low 32 bits of lane 0 come from NEG (whose low half is 0)
        if is_bv(x[0]):
            neg1 = lambda i, v: Concat(Extract(63, 32, v[0]), Extract(31, 0, NEG))
        else:  # clearing the low 32 bits of a two's complement value
            neg1 = lambda i, v: v[0] - v[0] % (1 << 32)
    else:
        neg1 = neg
    zero = lambda i, v: Z0
    m = [mx(x[i], y[i]) for i in range(4)]
    a = [sub(m[i], E[i]) for i in range(4)]
    b = [sub(Z0, add(Z[i], E[i])) for i in range(4)]
    if mutant == "perm":  # 0x94 = lanes (0, 1, 1, 2)
        t = [neg1(0, a), a[1], a[1], a[2]]
    else:
        t = shift(a, 1, neg1)
    a = [mx(a[i], add(t[i], b[i])) for i in range(4)]
    t = shift(b, 1, zero)
    b = [add(b[i], t[i]) for i in range(4)]
    t = shift(a, 2, neg)
    a = [mx(a[i], add(t[i], b[i])) for i in range(4)]
    t = shift(b, 2, zero)
    b = [add(b[i], t[i]) for i in range(4)]
    if mutant != "nocarry":
        a = [mx(a[i], add(carry, b[i])) for i in range(4)]
    Ip = [carry] + a[:3]
    V = []
    for i in range(4):
        v, _ = best3(x[i], y[i], sub(Ip[i], Z[i]))
        V.append(sub(v, E[i]))
    return a, Ip, V, [a[3]] * 4, ops


def scan_case(mutant=None, vlo=VLO, zlo=ZLO, elo=INT_MIN):
    """The block in exact integers, with the bounds, the serial recurrence it
    must equal, and every add/sub of the block and of 4 scalar steps."""
    x = [Int(f"x{i}") for i in range(4)]
    y = [Int(f"y{i}") for i in range(4)]
    Z = [Int(f"Z{i}") for i in range(4)]
    E = [Int(f"E{i}") for i in range(4)]
    carry = Int("carry")
    asm = [carry >= vlo, carry <= VHI]
    for i in range(4):
        asm += [x[i] >= vlo, x[i] <= VHI, y[i] >= vlo, y[i] <= VHI,
                Z[i] >= zlo, Z[i] <= ZHI, E[i] >= elo, E[i] <= INT_MAX]
    a, Ip, V, nc, ops = scan_lanes(x, y, Z, E, carry, mutant, const=IntVal)

    def mx(p, q):
        return If(p >= q, p, q)

    cl, I = [], carry
    for i in range(4):
        prev = I
        I = mx(mx(x[i], y[i]), prev - Z[i]) - E[i]
        cl += [a[i] == I, Ip[i] == prev, V[i] == I]
        # the scalar C code: BEST3's value, then minus E
        ops.append(("sub", prev, Z[i]))
        v, _ = best3(x[i], y[i], prev - Z[i])
        ops.append(("sub", v, E[i]))
        cl.append(v - E[i] == I)
    cl += [c == I for c in nc]
    return And(*cl), asm, ops


def in_range(ops):
    cl = []
    for op, p, q in ops:
        r = p + q if op == "add" else p - q
        cl += [r >= LLONG_MIN, r <= -LLONG_MIN - 1]
    return And(*cl)


def lemma_scan():
    claim, asm, _ = scan_case()
    return prove("L3 scan block (exact integers) == serial recurrence: I, I(col-1), stored, carry; "
                 "scalar step == recurrence", claim, asm)


def lemma_nowrap():
    _, asm, ops = scan_case()
    return prove(f"L4 all {len(ops)} adds/subs (block + 4 scalar steps) stay in 64-bit range",
                 in_range(ops), asm)


def lemma_scan_bounds():
    claim, asm, _ = scan_case(vlo=LLONG_MIN // 2 - (1 << 40))
    return prove("L3n with values allowed below the NEG fill, the block is wrong", claim, asm, expect=sat)


def self_test():
    ok = True
    for mut in ("perm", "negzero", "nocarry"):
        claim, asm, _ = scan_case(mutant=mut)
        ok &= prove(f"self-test mutant '{mut}' is refuted", claim, asm, expect=sat)
    return ok


def equiv():
    claim, asm, ops = scan_case(mutant="blend01", zlo=0, elo=0)
    ok = prove("equiv: step-1 NEG blend imm 0x01 equals 0x03 when Z >= 0 and E >= 0", claim, asm)
    ok &= prove("equiv: ... and its adds/subs stay in 64-bit range", in_range(ops), asm)
    claim, asm, _ = scan_case(mutant="blend01")
    ok &= prove("equiv: ... but not when Z + E may be negative", claim, asm, expect=sat)
    return ok


def main():
    args = sys.argv[1:]
    want = [a for a in args if not a.startswith("--")] or ["L0", "L1", "L2", "L3", "L4"]
    table = {"L0": lemma_intrinsics, "L1": lemma_best3, "L2": lemma_xy,
             "L3": lambda: lemma_scan() & lemma_scan_bounds(), "L4": lemma_nowrap}
    ok = True
    for name in want:
        ok &= table[name]()
    if "--self-test" in args:
        ok &= self_test()
    if "--equiv" in args:
        ok &= equiv()
    print("ALL OK" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


main()
