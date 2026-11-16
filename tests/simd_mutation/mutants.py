"""Mutants of the AVX2 kernels in mz_yama.c, for tests/simd_mutation/run_mutants.sh.

  python3 mutants.py list              names, one per line
  python3 mutants.py describe NAME     one-line description
  python3 mutants.py apply NAME FILE   rewrite FILE in place; exits 1 unless the
                                       old text occurs exactly once and the
                                       file changed

Each mutant is (name, description, old, new).  The old text must be unique
in mz_yama.c, so a mutant that stops matching after an edit fails loudly
instead of silently testing the unmutated code.
"""
import sys

MUTANTS = [
    ("perm1", "I scan step 1: permute a by 0x94 instead of 0x90",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x90), NEG, 0x03);",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x94), NEG, 0x03);"),
    ("perm2", "I scan step 2: permute a by 0x80 instead of 0x40",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x40), NEG, 0x0F);",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x80), NEG, 0x0F);"),
    ("permb", "I scan step 1: permute b by 0x94 instead of 0x90",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(b, 0x90), ZERO, 0x03);",
     "t = _mm256_blend_epi32(_mm256_permute4x64_epi64(b, 0x94), ZERO, 0x03);"),
    ("blend_ip", "I(col-1) blend takes the carry into lanes 0-1 (0x0F) instead of lane 0",
     "Ip = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x90), carry, 0x03);",
     "Ip = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x90), carry, 0x0F);"),
    ("blend01", "I scan step 1: NEG fills only the low 32 bits of lane 0 (0x01)",
     "_mm256_permute4x64_epi64(a, 0x90), NEG, 0x03);",
     "_mm256_permute4x64_epi64(a, 0x90), NEG, 0x01);"),
    ("negzero1", "I scan step 1: fill with ZERO instead of NEG",
     "_mm256_permute4x64_epi64(a, 0x90), NEG, 0x03);",
     "_mm256_permute4x64_epi64(a, 0x90), ZERO, 0x03);"),
    ("negzero2", "I scan step 2: fill with ZERO instead of NEG",
     "_mm256_permute4x64_epi64(a, 0x40), NEG, 0x0F);",
     "_mm256_permute4x64_epi64(a, 0x40), ZERO, 0x0F);"),
    ("nocarry", "I scan: the carried I(lo-1) is never applied",
     "        a = max_epi64(a, _mm256_add_epi64(carry, b));\n",
     "\n"),
    ("carrylane", "I scan: the new carry is lane 2 (0xAA) instead of lane 3",
     "carry = _mm256_permute4x64_epi64(a, 0xFF);",
     "carry = _mm256_permute4x64_epi64(a, 0xAA);"),
    ("noE", "I kernel: the stored value is not reduced by E",
     "        V = _mm256_sub_epi64(V, E);\n",
     "\n"),
    ("noEa", "I scan: a = m instead of m - E",
     "a = _mm256_sub_epi64(m, E);",
     "a = m;"),
    ("flagshift", "I kernel: flag_i packed at bit 5 instead of bit 4",
     "fb |= fi << 4;",
     "fb |= fi << 5;"),
    ("tie_di", "BEST3_V: D wins a D-I tie (BEST3 gives it to I)",
     "md = _mm256_and_si256(bx_, _mm256_cmpgt_epi64(y, z));",
     "md = _mm256_andnot_si256(_mm256_cmpgt_epi64(z, y), bx_);"),
    ("tie_c", "BEST3_V: C loses ties (BEST3 gives them to C)",
     "__m256i bx_ = _mm256_or_si256(_mm256_cmpgt_epi64(y, x), _mm256_cmpgt_epi64(z, x));",
     "__m256i bx_ = _mm256_xor_si256(_mm256_and_si256(_mm256_cmpgt_epi64(x, y), "
     "_mm256_cmpgt_epi64(x, z)), _mm256_set1_epi64x(-1));"),
    ("noDge", "C/D kernel: D is not reduced by the deleted column's extension",
     "Dv = _mm256_sub_epi64(Dv, Dge);",
     ";"),
    ("cshift", "C/D kernel: C reads pC[col] instead of pC[col-1]",
     "x = _mm256_loadu_si256((const __m256i *)(pC + col - 1));",
     "x = _mm256_loadu_si256((const __m256i *)(pC + col));"),
    ("p1start", "pass 1 vector start: LB[row-1/row-2] + 1 instead of + 2",
     "MAX(LB[row-1], LB[row-2]) + 2);",
     "MAX(LB[row-1], LB[row-2]) + 1);"),
    ("p1end", "pass 1 vector end: may include column N",
     "vhi = MIN(RB[row], N-1);",
     "vhi = MIN(RB[row], N);"),
    ("p2start", "pass 2 vector start: LB[row] + 1 instead of + 2",
     "vlo = LB[row] + 2;",
     "vlo = LB[row] + 1;"),
    ("p2end", "pass 2 vector end: rounds up past RB[row]",
     "vend = vlo + ((RB[row] - vlo + 1) & ~3);",
     "vend = vlo + ((RB[row] - vlo + 4) & ~3);"),
    # The check is the whole body of a brace-less for loop, so the mutant must
    # leave a statement there; an empty line made the next loop its body and
    # never terminated (the first run hung on it).
    ("fits_i_removed", "range check: Ix not checked",
     "                ok &= FITS32(r->Ix[k]);\n",
     "                ok &= 1;\n"),
    ("fits_cd_removed", "range check: Cx and Cz not checked",
     "ok &= FITS32(r->Cx[k]) & FITS32(r->Cz[k]);",
     "ok &= 1;"),
    ("fits_inverted", "range check inverted",
     "#define FITS32(v) ((v) >= INT_MIN && (v) <= INT_MAX)",
     "#define FITS32(v) (!((v) >= INT_MIN && (v) <= INT_MAX))"),
    # The value-bound check (YAMA_VBOUND etc.): these change no output (the bounds are never approached), so only
    # check.sh's "bound" cases, which read the shadow counters, can kill them.
    ("bound_ignored", "value-bound check: a failed check does not send the call to scalar",
     "simd = 0;\t// this call runs on the scalar path, exactly as with YAMA_SIMD=0",
     "(void)0;"),
    ("bound_v_removed", "value-bound check: (M+N)*Pmax not checked",
     " && (double)(M + N)*pmax <= YAMA_VBOUND", ""),
    ("bound_z_removed", "value-bound check: |Izc| not checked",
     "zmax <= YAMA_ZBOUND && ", ""),
    ("bound_e_removed", "value-bound check: K*L*|gap_extend| not checked",
     " && emax <= YAMA_EBOUND", ""),
]


def find(name):
    for m in MUTANTS:
        if m[0] == name:
            return m
    sys.exit(f"unknown mutant {name}")


def main():
    cmd = sys.argv[1]
    if cmd == "list":
        print("\n".join(m[0] for m in MUTANTS))
    elif cmd == "describe":
        print(find(sys.argv[2])[1])
    elif cmd == "apply":
        _, _, old, new = find(sys.argv[2])
        path = sys.argv[3]
        s = open(path).read()
        n = s.count(old)
        if n != 1:
            sys.exit(f"mutant {sys.argv[2]}: old text occurs {n} times in {path}, expected 1")
        t = s.replace(old, new)
        if t == s:
            sys.exit(f"mutant {sys.argv[2]}: file unchanged")
        open(path, "w").write(t)
    else:
        sys.exit(__doc__)


main()
