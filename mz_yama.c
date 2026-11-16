/* mz_yama.c -- dynamic programming procedure to align two alignments
*
* The algorithm uses sum-of-pairs substitution scores and the "quasi-natural
* gap costs" developed by:
*	Altschul, S. (1989) Gap costs for multiple sequence alignment.
*	J. theor. Biol. 138, 297-309.
* The graph model for aligning using affine gap scores, which employs an
* (M+1)-by-(N+1) grid (when the sequence lengths are M and N) with nodes called
* C, D and I at each grid point, is described by:
*	E. Myers and W. Miller (1989) Approximate matching of regular
*	expressions. Bull. Math. Biol. 51, 5-37.
* A slight modification to make the model appropriate for quasi-natural gap
* costs is explained by:
*	K.-M. Chao, R. Hardison and W. Miller (1994) Recent developments in
*	linear-space alignment methods: a survey. J. Comp. Biol. 1, 271-291.
* End-gaps (overhangs) are not charged a gap-open penalty.
*/

#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "mz_yama.h"

#define FLAG_C (0)
#define FLAG_I (1<<0)
#define FLAG_D (1<<1)
#define SELECT_CID (FLAG_I | FLAG_D | FLAG_C)
/* Path scores are 64-bit: a column scores the sum over all K*L row pairs, so
*  a long block of many rows passes INT_MAX (about 113,000 columns for 200
*  rows of identical sequence), which wrapped and silently misaligned the
*  rest of the block or sent the traceback off the grid.
*/
typedef long long score_t;
// MININT is a hugely negative number away from the underflow threshold
#define MININT (LLONG_MIN/4)

/* YAMA_ASSERT_BOUNDS (test build): abort if a value leaves the bounds the z3 proof of the I-pass scan assumes
*  (tests/simd_proof/scan_z3.py): stored C, D, I values in [VLO, VHI] = [-2^62 + 2^47, 2^62 - 2^47], |Izc| <= 2^44,
*  and eK = n*K*gap_extend within int.  It is checked where shadow mode already visits every value, so it implies
*  YAMA_SHADOW (each row is also compared cell by cell with the scalar code). */
#ifdef YAMA_ASSERT_BOUNDS
#ifndef YAMA_SHADOW
#define YAMA_SHADOW 1
#endif
#define YAMA_VHI ((1LL << 62) - (1LL << 47))
#define YAMA_VLO (-YAMA_VHI)
#define YAMA_ZMAX (1LL << 44)
#endif

/* Column counts.  Every per-cell sum in yama() runs over all K*L row pairs, but a
*  row of A enters a sum only through its symbol and whether it and the row
*  before it are dashes, and likewise a row of B.  So each sum is built from
*  counts of those classes: the gap sums from two 0/1 flags per row (4 classes),
*  the substitution sum from B's symbol counts.  The A side is folded into
*  per-row tables once, so a cell costs a few short dot products instead of
*  O(K*L) work.  Every term is an integer, so each total is exactly the value
*  the pair-by-pair loops produced, and the alignment is unchanged.
*/
#define CL(t,v) (((t)<<1)|(v))	// class of a row from two 0/1 dash flags

/* The best of the paths from the C, D and I nodes, with ties going to C, then D:
*  x if x >= y and x >= z, else y if y > z, else z. Written as selects rather than
*  branches, since which wins is unpredictable from cell to cell.
*/
#define BEST3(x, y, z, V, F) do { \
        int bx_ = ((x) >= (y)) & ((x) >= (z)), by_ = !bx_ & ((y) > (z)); \
        (V) = bx_ ? (x) : by_ ? (y) : (z); \
        (F) = bx_ ? FLAG_C : by_ ? FLAG_D : FLAG_I; \
    } while (0)
// a 4-class dot product: counts of B's row classes against one row table
#define DOT4(g,w) ((g)[0]*(w)[0] + (g)[1]*(w)[1] + (g)[2]*(w)[2] + (g)[3]*(w)[3])

typedef struct ROWTAB {
    score_t Ix[4], Iy[2], Cx[4], Cy[2], Cz[4], Dx[2], Dy, Dz[2];
}
row_tab;

/* Row-parallel SIMD for pass 1 (C and D, which read only the previous row).
*  YAMA_SIMD=0 at compile time keeps the scalar path only; otherwise, on x86-64
*  with GCC or Clang, an AVX2 kernel is compiled with a target attribute and
*  chosen at run time when the CPU has AVX2, so the default build still runs
*  on any x86-64.  The kernel does the band's interior, where every guard of
*  the scalar code is true; the edges stay scalar.  Its multiplies are exact
*  32x32->64 products, so a row whose coefficients do not all fit in 32 bits
*  is done by the scalar code.  Setting YAMA_SIMD=0 in the environment forces
*  the scalar path at run time (read once per process), for testing and as an
*  escape hatch; the output is the same either way.
*/
#ifndef YAMA_SIMD
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#define YAMA_SIMD 1
#else
#define YAMA_SIMD 0
#endif
#endif

#if YAMA_SIMD
#include <immintrin.h>
#include <string.h>

// SPREAD[m] puts bit j of the 4-bit mask m into bit 0 of byte j
static const unsigned SPREAD[16] = {
    0x00000000, 0x00000001, 0x00000100, 0x00000101,
    0x00010000, 0x00010001, 0x00010100, 0x00010101,
    0x01000000, 0x01000001, 0x01000100, 0x01000101,
    0x01010000, 0x01010001, 0x01010100, 0x01010101
};

#ifndef FITS32	// overridable, to test the scalar fallback
#define FITS32(v) ((v) >= INT_MIN && (v) <= INT_MAX)
#endif

/* The value bounds the z3 proof of the I scan assumes (tests/simd_proof/scan_z3.py): every C, D, I value and
*  every candidate the kernels form within [-2^62 + 2^47, 2^62 - 2^47], |Izc| <= 2^44, and eK within int.  yama()
*  checks a sufficient condition once per call, before the first row, in O(1), and runs the whole call on the scalar path
*  if it fails.  The limits are overridable, so a test build can lower one and see the fallback taken
*  (tests/yama_fuzz/check.sh, "bound"). */
#ifndef YAMA_VBOUND	// (M+N) times a bound on the change one DP edge can make; see yama()
#define YAMA_VBOUND ((double)(1LL << 60))
#endif
#ifndef YAMA_ZBOUND	// K*L*max|gop|, which bounds every |Izc[col]|
#define YAMA_ZBOUND ((double)(1LL << 44))
#endif
#ifndef YAMA_EBOUND	// K*L*|gap_extend|, which bounds every eK and every D-edge extension
#define YAMA_EBOUND ((double)INT_MAX)
#endif
static double dabs(double v) {
    return v < 0 ? -v : v;
}

/* BEST3 on 4 lanes: V is the winner, mi and md flag lanes won by z (I) and
*  y (D), so the 2-bit flag is FLAG_I*mi | FLAG_D*md, exactly as BEST3.
*/
#define BEST3_V(x, y, z, V, mi, md) do { \
        __m256i bx_ = _mm256_or_si256(_mm256_cmpgt_epi64(y, x), _mm256_cmpgt_epi64(z, x)); \
        /* bx_ is NOT (x >= y && x >= z) */ \
        md = _mm256_and_si256(bx_, _mm256_cmpgt_epi64(y, z)); \
        mi = _mm256_andnot_si256(md, bx_); \
        V = _mm256_blendv_epi8(_mm256_blendv_epi8(x, y, md), z, mi); \
    } while (0)

/* Columns lo <= col < hi (hi - lo a multiple of 4) of pass 1 for one row in
*  the band's interior: row > 1, col >= 2, col < N, col >= LB[row-1]+2 and
*  col >= LB[row-2]+2.  cs holds per-column int64 copies of gb[.][0..3],
*  vb[.][0..1] and hist[.][0..nsym-1], each stride long.
*/
__attribute__((target("avx2")))
static void cd_avx2(int lo, int hi, const score_t *pC, const score_t *pD,
                    const score_t *pI, score_t *cC, score_t *cD, uchar *fcd,
                    const long long *cs, size_t stride, int nsym,
                    const row_tab *r, const score_t *w, score_t dge) {
    const long long *g0 = cs, *g1 = cs + stride, *g2 = cs + 2*stride,
                     *g3 = cs + 3*stride, *v0 = cs + 4*stride, *v1 = cs + 5*stride,
                      *h = cs + 6*stride;
    __m256i Cx0 = _mm256_set1_epi64x(r->Cx[0]), Cx1 = _mm256_set1_epi64x(r->Cx[1]),
            Cx2 = _mm256_set1_epi64x(r->Cx[2]), Cx3 = _mm256_set1_epi64x(r->Cx[3]),
            Cz0 = _mm256_set1_epi64x(r->Cz[0]), Cz1 = _mm256_set1_epi64x(r->Cz[1]),
            Cz2 = _mm256_set1_epi64x(r->Cz[2]), Cz3 = _mm256_set1_epi64x(r->Cz[3]),
            Cy0 = _mm256_set1_epi64x(r->Cy[0]), Cy1 = _mm256_set1_epi64x(r->Cy[1]),
            Dx0 = _mm256_set1_epi64x(r->Dx[0]), Dx1 = _mm256_set1_epi64x(r->Dx[1]),
            Dz0 = _mm256_set1_epi64x(r->Dz[0]), Dz1 = _mm256_set1_epi64x(r->Dz[1]),
            Dy = _mm256_set1_epi64x(r->Dy), Dge = _mm256_set1_epi64x(dge);
    __m256i W[16];
    int col, k, nw = nsym < 16 ? nsym : 16;

    for (k = 0; k < nw; ++k)
        W[k] = _mm256_set1_epi64x(w[k]);
    for (col = lo; col < hi; col += 4) {
        __m256i G0 = _mm256_loadu_si256((const __m256i *)(g0 + col)),
                G1 = _mm256_loadu_si256((const __m256i *)(g1 + col)),
                G2 = _mm256_loadu_si256((const __m256i *)(g2 + col)),
                G3 = _mm256_loadu_si256((const __m256i *)(g3 + col)),
                V0 = _mm256_loadu_si256((const __m256i *)(v0 + col)),
                V1 = _mm256_loadu_si256((const __m256i *)(v1 + col));
        __m256i x, y, z, Cv, Dv, ic, dc, id, dd, t;
        unsigned b;

        // C from (row-1,col-1)
        x = _mm256_loadu_si256((const __m256i *)(pC + col - 1));
        y = _mm256_loadu_si256((const __m256i *)(pD + col - 1));
        z = _mm256_loadu_si256((const __m256i *)(pI + col - 1));
        t = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epi32(G0, Cx0), _mm256_mul_epi32(G1, Cx1)),
                             _mm256_add_epi64(_mm256_mul_epi32(G2, Cx2), _mm256_mul_epi32(G3, Cx3)));
        x = _mm256_sub_epi64(x, t);
        t = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epi32(G0, Cz0), _mm256_mul_epi32(G1, Cz1)),
                             _mm256_add_epi64(_mm256_mul_epi32(G2, Cz2), _mm256_mul_epi32(G3, Cz3)));
        z = _mm256_sub_epi64(z, t);
        t = _mm256_add_epi64(_mm256_mul_epi32(V0, Cy0), _mm256_mul_epi32(V1, Cy1));
        y = _mm256_sub_epi64(y, t);
        BEST3_V(x, y, z, Cv, ic, dc);
        for (k = 0; k < nw; ++k)
            Cv = _mm256_add_epi64(Cv, _mm256_mul_epi32(
                     _mm256_loadu_si256((const __m256i *)(h + k*stride + col)), W[k]));
        for ( ; k < nsym; ++k)
            Cv = _mm256_add_epi64(Cv, _mm256_mul_epi32(
                     _mm256_loadu_si256((const __m256i *)(h + k*stride + col)),
                     _mm256_set1_epi64x(w[k])));

        // D from (row-1,col)
        x = _mm256_loadu_si256((const __m256i *)(pC + col));
        y = _mm256_loadu_si256((const __m256i *)(pD + col));
        z = _mm256_loadu_si256((const __m256i *)(pI + col));
        x = _mm256_sub_epi64(x, _mm256_add_epi64(_mm256_mul_epi32(V0, Dx0), _mm256_mul_epi32(V1, Dx1)));
        y = _mm256_sub_epi64(y, Dy);
        z = _mm256_sub_epi64(z, _mm256_add_epi64(_mm256_mul_epi32(V0, Dz0), _mm256_mul_epi32(V1, Dz1)));
        BEST3_V(x, y, z, Dv, id, dd);
        Dv = _mm256_sub_epi64(Dv, Dge);

        _mm256_storeu_si256((__m256i *)(cC + col), Cv);
        _mm256_storeu_si256((__m256i *)(cD + col), Dv);
        // fcd = flag_c | flag_d<<2, with FLAG_I = 1 and FLAG_D = 2
        b = SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(ic))]
            | SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(dc))] << 1
            | SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(id))] << 2
            | SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(dd))] << 3;
        memcpy(fcd + col, &b, 4);	// little-endian: byte j is column col+j
    }
}

/* Pass 2 (I) for columns lo <= col < hi (hi - lo a multiple of 4) of a row
*  below M, in the interior col >= LB[row]+2, where both of its guards hold.
*  I is a running maximum, so it is found by a prefix scan in the (max,+)
*  semiring, four columns at a time.  BEST3's value is max(x, y, z) whatever
*  its tie order, so with m = max(x, y), e = vb[col][0]*K*gap_extend and
*  Z = Izc[col]:
*      I(col) = max(m, I(col-1) - Z) - e = max(a, I(col-1) + b),
*  a = m - e, b = -(Z + e).  (a2,b2) o (a1,b1) = (max(a2, a1 + b2), b1 + b2)
*  is associative, so a 2-step Hillis-Steele scan over the 4 lanes and the
*  carried I(lo-1) give I for each lane exactly (integers, far from overflow:
*  MININT is LLONG_MIN/4).  With I(col-1) known in every lane, BEST3 then runs
*  as in the scalar code for the flags, and its value less e is I(col).
*/
#ifndef YAMA_SIMD_I
#define YAMA_SIMD_I 1
#endif
#if YAMA_SIMD_I
__attribute__((target("avx2")))
static inline __m256i max_epi64(__m256i a, __m256i b) {
    return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(b, a));
}

__attribute__((target("avx2")))
static score_t i_avx2(int lo, int hi, const score_t *cC, const score_t *cD,
                      score_t *cI, const uchar *fcd, uchar *tbp,
                      const long long *cs, size_t stride, const score_t *Izc,
                      const long long *eK, const row_tab *r, score_t Iin) {
    const long long *g0 = cs, *g1 = cs + stride, *g2 = cs + 2*stride,
                     *g3 = cs + 3*stride, *v0 = cs + 4*stride, *v1 = cs + 5*stride;
    const __m256i Ix0 = _mm256_set1_epi64x(r->Ix[0]), Ix1 = _mm256_set1_epi64x(r->Ix[1]),
                  Ix2 = _mm256_set1_epi64x(r->Ix[2]), Ix3 = _mm256_set1_epi64x(r->Ix[3]),
                  Iy0 = _mm256_set1_epi64x(r->Iy[0]), Iy1 = _mm256_set1_epi64x(r->Iy[1]),
                  NEG = _mm256_set1_epi64x(LLONG_MIN/2), ZERO = _mm256_setzero_si256();
    __m256i carry = _mm256_set1_epi64x(Iin);
    int col;

    for (col = lo; col < hi; col += 4) {
        __m256i G0 = _mm256_loadu_si256((const __m256i *)(g0 + col)),
                G1 = _mm256_loadu_si256((const __m256i *)(g1 + col)),
                G2 = _mm256_loadu_si256((const __m256i *)(g2 + col)),
                G3 = _mm256_loadu_si256((const __m256i *)(g3 + col)),
                V0 = _mm256_loadu_si256((const __m256i *)(v0 + col)),
                V1 = _mm256_loadu_si256((const __m256i *)(v1 + col)),
                Z = _mm256_loadu_si256((const __m256i *)(Izc + col)),
                E = _mm256_loadu_si256((const __m256i *)(eK + col));
        __m256i x, y, z, m, a, b, t, Ip, V, mi, md;
        unsigned fb, fi;

        x = _mm256_loadu_si256((const __m256i *)(cC + col - 1));
        y = _mm256_loadu_si256((const __m256i *)(cD + col - 1));
        t = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epi32(G0, Ix0), _mm256_mul_epi32(G1, Ix1)),
                             _mm256_add_epi64(_mm256_mul_epi32(G2, Ix2), _mm256_mul_epi32(G3, Ix3)));
        x = _mm256_sub_epi64(x, t);
        y = _mm256_sub_epi64(y, _mm256_add_epi64(_mm256_mul_epi32(V0, Iy0), _mm256_mul_epi32(V1, Iy1)));
        m = max_epi64(x, y);
        a = _mm256_sub_epi64(m, E);
        b = _mm256_sub_epi64(ZERO, _mm256_add_epi64(Z, E));
        // step 1: combine with the lane to the left
        t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x90), NEG, 0x03);
        a = max_epi64(a, _mm256_add_epi64(t, b));
        t = _mm256_blend_epi32(_mm256_permute4x64_epi64(b, 0x90), ZERO, 0x03);
        b = _mm256_add_epi64(b, t);
        // step 2: combine with the lane two to the left
        t = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x40), NEG, 0x0F);
        a = max_epi64(a, _mm256_add_epi64(t, b));
        t = _mm256_blend_epi32(_mm256_permute4x64_epi64(b, 0x40), ZERO, 0x0F);
        b = _mm256_add_epi64(b, t);
        // apply the carried I(lo-1): a is now I for each lane
        a = max_epi64(a, _mm256_add_epi64(carry, b));
        // I(col-1) in each lane, then BEST3 exactly as the scalar code
        Ip = _mm256_blend_epi32(_mm256_permute4x64_epi64(a, 0x90), carry, 0x03);
        carry = _mm256_permute4x64_epi64(a, 0xFF);
        z = _mm256_sub_epi64(Ip, Z);
        BEST3_V(x, y, z, V, mi, md);
        V = _mm256_sub_epi64(V, E);
        _mm256_storeu_si256((__m256i *)(cI + col), V);
        fi = SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(mi))]
             | SPREAD[_mm256_movemask_pd(_mm256_castsi256_pd(md))] << 1;
        memcpy(&fb, fcd + col, 4);
        fb |= fi << 4;
        memcpy(tbp + (col - lo), &fb, 4);
    }
    return cI[hi-1];
}
#endif
#endif

#ifdef YAMA_SHADOW
/* Shadow mode (test instrument): every row is recomputed by a copy of the
*  scalar code into separate buffers and compared cell by cell with what the
*  production path (vector kernels plus scalar edges) stored: C, D and I
*  values and every traceback flag.  The first mismatch aborts.  Counters,
*  printed at exit, show how many cells were compared and how many rows went
*  through each kernel, so a run that compared nothing cannot pass silently.
*/
#include <stdio.h>
static unsigned long long sh_calls, sh_rows, sh_cells, sh_p1_vrows, sh_p1_vcells,
       sh_p2_vrows, sh_p2_vcells, sh_unfit_rows, sh_bound_fallbacks;
static score_t sh_vmin = LLONG_MAX, sh_vmax = LLONG_MIN, sh_vmin_nonsentinel = LLONG_MAX,
               sh_zmax, sh_emax;
static void sh_report(void) {
    /* YAMA_SHADOW_LOG names a file to append the counters to as well, for
    *  drivers that discard stderr (tests/run.sh, roast's children). */
    const char *lp = getenv("YAMA_SHADOW_LOG");
    FILE *f = lp ? fopen(lp, "a") : NULL;
    int i;
    for (i = 0; i < 2; ++i) {
        FILE *o = i ? f : stderr;
        if (!o)
            continue;
        fprintf(o, "YAMA_SHADOW calls=%llu rows=%llu cells=%llu p1_vec_rows=%llu "
            "p1_vec_cells=%llu p2_vec_rows=%llu p2_vec_cells=%llu unfit_rows=%llu "
            "vmin=%lld vmin_above_MININT=%lld vmax=%lld max_abs_Izc=%lld max_abs_eK=%lld "
            "bound_fallbacks=%llu\n",
            sh_calls, sh_rows, sh_cells, sh_p1_vrows, sh_p1_vcells, sh_p2_vrows,
            sh_p2_vcells, sh_unfit_rows, sh_vmin, sh_vmin_nonsentinel, sh_vmax,
            sh_zmax, sh_emax, sh_bound_fallbacks);
    }
    if (f)
        fclose(f);
}
static void sh_report(void);
static void sh_val(score_t v) {
#ifdef YAMA_ASSERT_BOUNDS
    if (v < YAMA_VLO || v > YAMA_VHI) {
        fprintf(stderr, "YAMA_ASSERT_BOUNDS value %lld outside [%lld, %lld]\n", (long long)v,
                (long long)YAMA_VLO, (long long)YAMA_VHI);
        sh_report();
        abort();
    }
#endif
    if (v < sh_vmin) sh_vmin = v;
    if (v > sh_vmax) sh_vmax = v;
    if (v > MININT && v < sh_vmin_nonsentinel) sh_vmin_nonsentinel = v;
}
static void sh_fail(const char *field, int row, int col, long long got, long long want,
                    int lb, int rb, int M, int N, int K, int L, int p1lo, int p1end,
                    int p2lo, int p2end) {
    fprintf(stderr, "YAMA_SHADOW MISMATCH field=%s row=%d col=%d production=%lld scalar=%lld "
            "band=[%d,%d] M=%d N=%d K=%d L=%d p1_vec=[%d,%d) p2_vec=[%d,%d)\n",
            field, row, col, got, want, lb, rb, M, N, K, L, p1lo, p1end, p2lo, p2end);
    if (getenv("YAMA_SHADOW_LOG")) {
        FILE *f = fopen(getenv("YAMA_SHADOW_LOG"), "a");
        if (f) {
            fprintf(f, "YAMA_SHADOW MISMATCH field=%s row=%d col=%d production=%lld scalar=%lld\n",
                    field, row, col, got, want);
            fclose(f);
        }
    }
    sh_report();
    abort();
}
#define SH_CMP(field, got, want) do { \
        if ((long long)(got) != (long long)(want)) \
            sh_fail(field, row, col, (long long)(got), (long long)(want), LB[row], RB[row], \
                    M, N, K, L, sh_p1lo, sh_p1end, vlo, vend); \
    } while (0)
#endif

// populate a column from a column of height K and a column of height L
static void new_col(uchar *from1, int K, uchar *from2, int L, uchar *to_col) {
    //void new_col(uchar *from1, int K, uchar *from2, int L, uchar *to_col) {
    int i;

    for (i = 0; i < K; ++i)
        to_col[i] = from1[i];
    for (i = 0; i < L; ++i)
        to_col[i+K] = from2[i];
}

// Entry point: (documented in mz_yama.h)
void yama(uchar **A, int K, int M, uchar **B, int L, int N, int *LB, int *RB,
          uchar ***OAL, int *OM) {
    int i, j, k, m, m_new, n, row, col, s, t, u, v, nedit, nsym, c;
    score_t C, D, I, x, y, z, Iz[4];
    int (*gb)[4], (*vb)[2], (*cnt)[4], (*ua)[2], *hist, sym[256], sym_ix[256];
    score_t *wA, *Izc;
    row_tab *rt;
    size_t tback_size;    // sum of band widths; can exceed INT_MAX for long insertions
    uchar **tback_row, *tback, *tbp, st, flag_d, flag_i, flag_c, node,
    *dashes, *script, **al_new, *fcd;
    score_t *dpbuf, *pC, *pD, *pI, *cC, *cD, *cI, *sw;
    int vlo, vend;
#if YAMA_SIMD
    int simd, vhi, *rowok = NULL;
    long long *cs = NULL;
    size_t stride = 0;
    double zmax, emax, pmax;	// the value-bound check
#endif
#if YAMA_SIMD && YAMA_SIMD_I
    long long *eK = NULL;
#endif
#ifdef YAMA_SHADOW
    score_t *sC, *sD, sx, sy, sz, sV;
    uchar *sfcd;
    int sf, sg, sk, sh_p1lo, sh_p1end;
    static int sh_registered;
#endif

#ifdef YAMA_DUMP_CASES
    /* Test build: append this call's inputs, in tests/yama_fuzz/main.c's case format, to the file named by
    *  $YAMA_DUMP_CASES, so real calls (from multiz on real data) become a fixed corpus for the fuzz check. */
    {
        const char *dp = getenv("YAMA_DUMP_CASES");
        FILE *df = dp ? fopen(dp, "a") : NULL;
        if (df) {
            fprintf(df, "%d %d\n", K, M);
            for (i = 1; i <= M; i++) fprintf(df, "%.*s\n", K, A[i]);
            fprintf(df, "%d %d\n", L, N);
            for (i = 1; i <= N; i++) fprintf(df, "%.*s\n", L, B[i]);
            for (i = 0; i <= M; i++) fprintf(df, "%d%c", LB[i], i == M ? '\n' : ' ');
            for (i = 0; i <= M; i++) fprintf(df, "%d%c", RB[i], i == M ? '\n' : ' ');
            fclose(df);
        }
    }
#endif

    if (LB[0] != 0 || RB[M] != N)
        fatalf("LB and RB not terminated properly: %d %d %d", LB[0], RB[M], N);
    for (tback_size = row = 0; row <= M; row++) {
        j = RB[row] - LB[row];
        // temporary sanity check:
        if (j < MIN(N,10)) {
            fatalf("RB[%d] - LB[%d] < %d, %d %d %d", row, row, MIN(N,10), RB[row], LB[row], N);
        }
        tback_size += (size_t)(j+1);
        if (row > 0 && LB[row] < LB[row-1])
            fatal("LB not monotonic");
        if (row > 0 && RB[row] < RB[row-1])
            fatal("RB not monotonic");
    }

    dashes = ckalloc(MAX(K,L)*sizeof(uchar));
    for (i = 0; i < MAX(K,L); ++i)
        dashes[i] = '-';

    tback = ckalloc(tback_size * sizeof(uchar));
    tback_row = ckalloc((M+1)*sizeof(uchar *));
    tbp = tback_row[0] = tback;
    *(tbp++) = 0;	/* unused */

    /* Two DP rows, previous (p) and current (c), stored as separate C, D and
    *  I arrays.  Columns of p beyond the previous row's band are never written
    *  by any row (RB is monotonic), so they hold MININT as the single in-place
    *  row did; columns below the band are never read.
    */
    pC = dpbuf = ckalloc(6*(size_t)(N+1)*sizeof(score_t));
    pD = pC + (N+1);
    pI = pD + (N+1);
    cC = pI + (N+1);
    cD = cC + (N+1);
    cI = cD + (N+1);
    fcd = ckalloc((N+1)*sizeof(uchar));
#ifdef YAMA_SHADOW
    if (!sh_registered) {
        atexit(sh_report);
        sh_registered = 1;
    }
    ++sh_calls;
    sC = ckalloc(2*(size_t)(N+1)*sizeof(score_t));
    sD = sC + (N+1);
    sfcd = ckalloc((N+1)*sizeof(uchar));
#endif
    for (col = 0; col <= N; ++col)
        cC[col] = cD[col] = cI[col] = MININT;
    pC[0] = pD[0] = pI[0] = 0;
    flag_i = FLAG_I;
    for (col = 1; col <= RB[0]; ++col) {
        pC[col] = pD[col] = MININT;
        for (j = n = 0; j < L; ++j)
            if (B[col][j] != '-')
                ++n;
        pI[col] = pI[col-1] - n*K*gap_extend;
        *(tbp++) = (flag_i << 4);
    }
    for ( ; col <= N; ++col)
        pC[col] = pD[col] = pI[col] = MININT;

    /* B side: for column col, gb counts rows by (B[col-1] is '-', B[col] is
    *  '-'), vb by B[col] alone, and hist counts each symbol.
    */
    /* Symbols whose substitution-score columns are identical (upper and lower
    *  case, and every non-ACGT letter, under the shipped matrices) share one
    *  class: their terms in the per-cell sum have equal weights, so summing
    *  their counts first gives the same total with fewer terms.
    */
    for (c = 0; c < 256; ++c)
        sym_ix[c] = -1;
    for (nsym = 0, col = 1; col <= N; ++col)
        for (j = 0; j < L; ++j)
            if (sym_ix[B[col][j]] < 0) {
                c = B[col][j];
                for (k = 0; k < nsym; ++k) {
                    for (i = 0; i < SS_NCHARS && SS(i, c) == SS(i, sym[k]); ++i)
                        ;
                    if (i == SS_NCHARS)
                        break;
                }
                sym_ix[c] = k;
                if (k == nsym)
                    sym[nsym++] = c;
            }
    gb = ckalloc((N+1)*sizeof(*gb));
    vb = ckalloc((N+1)*sizeof(*vb));
    hist = ckalloc((size_t)(N+1)*(nsym ? nsym : 1)*sizeof(int));
    for (col = 1; col <= N; ++col) {
        gb[col][0] = gb[col][1] = gb[col][2] = gb[col][3] = 0;
        vb[col][0] = vb[col][1] = 0;
        for (k = 0; k < nsym; ++k)
            hist[(size_t)col*nsym + k] = 0;
        for (j = 0; j < L; ++j) {
            t = (col > 1 && B[col-1][j] == '-');
            v = (B[col][j] == '-');
            gb[col][CL(t,v)]++;
            vb[col][v]++;
            hist[(size_t)col*nsym + sym_ix[B[col][j]]]++;
        }
    }

    /* A side: for row, cnt counts rows by (A[row-1] is '-', A[row] is '-'),
    *  ua by A[row] alone.  rt holds each gap sum's A half against every B
    *  class, and wA[row][k] is the sum over A[row] of ss[.][sym[k]].
    */
    cnt = ckalloc((M+1)*sizeof(*cnt));
    ua = ckalloc((M+1)*sizeof(*ua));
    rt = ckalloc((M+1)*sizeof(row_tab));
    wA = ckalloc((size_t)(M+1)*(nsym ? nsym : 1)*sizeof(score_t));
    for (t = 0; t < 2; ++t)
        for (v = 0; v < 2; ++v)
            Iz[CL(t,v)] = (score_t)K*GAP(1,t,1,v);
    for (row = 1; row <= M; ++row) {
        cnt[row][0] = cnt[row][1] = cnt[row][2] = cnt[row][3] = 0;
        ua[row][0] = ua[row][1] = 0;
        for (k = 0; k < nsym; ++k)
            wA[(size_t)row*nsym + k] = 0;
        for (i = 0; i < K; ++i) {
            s = (row > 1 && A[row-1][i] == '-');
            u = (A[row][i] == '-');
            cnt[row][CL(s,u)]++;
            ua[row][u]++;
            for (k = 0; k < nsym; ++k)
                wA[(size_t)row*nsym + k] += SS(A[row][i],sym[k]);
        }
        rt[row].Dy = 0;
        for (t = 0; t < 2; ++t) {
            rt[row].Dx[t] = rt[row].Dz[t] = 0;
            rt[row].Iy[t] = rt[row].Cy[t] = 0;
            for (v = 0; v < 2; ++v)
                rt[row].Ix[CL(t,v)] = rt[row].Cx[CL(t,v)] = rt[row].Cz[CL(t,v)] = 0;
        }
        for (s = 0; s < 2; ++s)
            for (u = 0; u < 2; ++u) {
                n = cnt[row][CL(s,u)];
                rt[row].Dy += (score_t)n*L*GAP(s,1,u,1);
                for (t = 0; t < 2; ++t) {
                    rt[row].Dx[t] += (score_t)n*GAP(s,t,u,1);
                    rt[row].Cy[t] += (score_t)n*GAP(s,1,u,t);	// t plays v here
                    for (v = 0; v < 2; ++v)
                        rt[row].Cx[CL(t,v)] += (score_t)n*GAP(s,t,u,v);
                }
            }
        for (s = 0; s < 2; ++s)	// I: s is A[row]'s own dash flag, u is 1
            for (t = 0; t < 2; ++t) {
                rt[row].Iy[t] += (score_t)ua[row][s]*GAP(s,1,1,t);	// t plays v
                for (v = 0; v < 2; ++v)
                    rt[row].Ix[CL(t,v)] += (score_t)ua[row][s]*GAP(s,t,1,v);
            }
        for (u = 0; u < 2; ++u)
            for (t = 0; t < 2; ++t) {
                rt[row].Dz[t] += (score_t)ua[row][u]*GAP(1,t,u,1);
                for (v = 0; v < 2; ++v)
                    rt[row].Cz[CL(t,v)] += (score_t)ua[row][u]*GAP(1,t,u,v);
            }
    }

    // The I node's own gap sum depends on the column alone.
    Izc = ckalloc((N+1)*sizeof(score_t));
    for (col = 1; col <= N; ++col)
        Izc[col] = DOT4(gb[col], Iz);
#ifdef YAMA_SHADOW
    for (col = 1; col <= N; ++col) {
        sV = Izc[col] < 0 ? -Izc[col] : Izc[col];
        if (sV > sh_zmax) sh_zmax = sV;
        n = vb[col][0];
        sV = (score_t)n*K*gap_extend;   /* in 64 bits, so an int overflow is seen, not committed */
        if (sV < 0) sV = -sV;
        if (sV > sh_emax) sh_emax = sV;
#ifdef YAMA_ASSERT_BOUNDS
        if (sh_zmax > YAMA_ZMAX || sh_emax > INT_MAX) {
            fprintf(stderr, "YAMA_ASSERT_BOUNDS |Izc| %lld (max %lld) or |eK| %lld (max %d) out of bounds\n",
                    (long long)sh_zmax, (long long)YAMA_ZMAX, (long long)sh_emax, INT_MAX);
            sh_report();
            abort();
        }
#endif
    }
#endif

#if YAMA_SIMD
    {   /* run-time choice: the CPU has AVX2 and YAMA_SIMD=0 is not set */
        static int allowed = -1;
        if (allowed < 0) {
            const char *e = getenv("YAMA_SIMD");
            allowed = !(e && e[0] == '0' && e[1] == '\0');
        }
        simd = allowed && __builtin_cpu_supports("avx2");
    }
    if (simd) {
        stride = (size_t)N + 1;
        cs = ckalloc((6 + (size_t)nsym)*stride*sizeof(long long));
        for (col = 1; col <= N; ++col) {
            for (k = 0; k < 4; ++k)
                cs[k*stride + col] = gb[col][k];
            cs[4*stride + col] = vb[col][0];
            cs[5*stride + col] = vb[col][1];
            for (k = 0; k < nsym; ++k)
                cs[(6 + k)*stride + col] = hist[(size_t)col*nsym + k];
        }
        rowok = ckalloc((M+1)*sizeof(int));
        /* The value-bound check (see YAMA_VBOUND), in O(1) per call.  Every DP value is reached from (0,0), or from a
        *  MININT sentinel, along at most M+N edges, and one edge changes a value by at most
        *      P = K*L*(gmax + smax) + K*L*gmax + K*L*gmax + K*L*|gap_extend|,
        *  with gmax = max|gop| and smax = max|ss|.  Each row coefficient (Cx, Cy, Cz, Dx, Dz, Ix, Iy) is a sum over A's
        *  K rows of counts times gop entries, so at most K*gmax; Dy and every Izc are at most K*L*gmax; wA is a sum of
        *  K ss entries; and the counts of B's rows that multiply them sum to L.  So the C, D or I term of an edge is at
        *  most K*L*gmax, its substitution sum K*L*smax, Dy or Izc K*L*gmax, and its extension K*L*|gap_extend|.
        *  Hence every value and candidate lies in [MININT - (M+N)*P, (M+N)*P].  With (M+N)*P <= 2^60 and
        *  MININT = -2^61 that is inside [-2^62 + 2^47, 2^62 - 2^47], and |Izc| <= zmax = K*L*gmax <= 2^44 and eK within
        *  int follow too.  The sums are in double, whose relative error here (a few roundings of 2^-53) cannot carry a
        *  bound of 2^60 past 2^61 - 2^47, and which cannot overflow.  Under the shipped tables (gap open at most 600,
        *  |ss| at most 157, gap_extend at most 50) and K = L = 216, P < 2^27, and M + N < 2^32, so the check cannot
        *  fail on such input; it stops a larger table or a far deeper input from leaving the proved range unnoticed.
        *  (A per-row form, from each row's actual coefficients, cost 3.5-4.2% more instructions on whole runs.) */
        {
            static int **ss_seen;	// ss tables are built once (init_scores70/85) and never changed, so cache by table
            static double smax;
            double gmax = 0;
            if (ss != ss_seen) {
                smax = 0;
                for (i = 0; i < SS_NCHARS; ++i)
                    for (j = 0; j < SS_NCHARS; ++j)
                        smax = MAX(smax, dabs((double)SS(i, j)));
                ss_seen = ss;
            }
            for (k = 0; k < 16; ++k)	// gop can be changed between calls (the fuzz driver does), so read it each time
                gmax = MAX(gmax, dabs((double)gop[k]));
            zmax = (double)K*L*gmax;
            emax = dabs((double)K*L*gap_extend);
            pmax = (double)K*L*(gmax + smax) + 2*zmax + emax;
        }
        /* The FITS32 checks keep _mm256_mul_epi32, which multiplies the low 32 bits only, exact.  Under the shipped
        *  gap table every nonzero gap-open entry is equal, so the Cx, Cz and Ix checks are implied by the Cy, Dx, Dz
        *  and Iy checks: tests/simd_mutation/real_table_checks.sh removes them at every magnitude of a real-shaped
        *  table and the output never changes.  They are kept as a defence anyway, at a few comparisons per row:
        *  they guard any other gap table (a configurable one, or a test's) and any input that breaks the implication,
        *  and check.sh's "fits32 idx" cases, which set one gap-open entry at a time, fail if any one is removed. */
        for (row = 1; row <= M; ++row) {
            row_tab *r = rt + row;
            int ok = 1;
            for (k = 0; k < 4; ++k)
                ok &= FITS32(r->Cx[k]) & FITS32(r->Cz[k]);
            for (k = 0; k < 2; ++k)
                ok &= FITS32(r->Cy[k]) & FITS32(r->Dx[k]) & FITS32(r->Dz[k]);
            for (k = 0; k < nsym; ++k)
                ok &= FITS32(wA[(size_t)row*nsym + k]);
#if YAMA_SIMD_I
            for (k = 0; k < 4; ++k)
                ok &= FITS32(r->Ix[k]);
            for (k = 0; k < 2; ++k)
                ok &= FITS32(r->Iy[k]);
#endif
            rowok[row] = ok;
        }
        if (!(zmax <= YAMA_ZBOUND && emax <= YAMA_EBOUND && (double)(M + N)*pmax <= YAMA_VBOUND)) {
            simd = 0;	// this call runs on the scalar path, exactly as with YAMA_SIMD=0
#ifdef YAMA_SHADOW
            ++sh_bound_fallbacks;
#endif
        }
    }
    if (simd) {
#if YAMA_SIMD_I
        eK = ckalloc(stride*sizeof(long long));
        for (col = 1; col <= N; ++col) {
            n = vb[col][0];
            eK[col] = n*K*gap_extend;	// the int product, as the scalar code
        }
#endif
    }
#endif

    C = D = I = MININT;
    for (row = 1; row <= M; ++row) {
        tback_row[row] = tbp - LB[row];

        /* Pass 1: C and D for every column of the band.  C(row,col) reads only
        *  (row-1,col-1) and D(row,col) only (row-1,col), so the columns are
        *  independent.
        */
        vlo = vend = -1;	// no vector segment
#if !YAMA_SIMD
        (void)vlo; (void)vend;	// read only by the vector and shadow code; keeps scalar builds warning-free
#endif
#if YAMA_SIMD
        if (simd && row > 1 && rowok[row]) {
            // the interior, where every guard below is true
            vlo = MAX(MAX(LB[row], 2), MAX(LB[row-1], LB[row-2]) + 2);
            vhi = MIN(RB[row], N-1);
            vend = vlo + ((vhi - vlo + 1) & ~3);
            if (vend <= vlo)
                vlo = vend = -1;
        }
#endif
#ifdef YAMA_SHADOW
        sh_p1lo = vlo;
        sh_p1end = vend;
        ++sh_rows;
#if YAMA_SIMD
        if (simd && !rowok[row])
            ++sh_unfit_rows;
#endif
        if (vlo >= 0) {
            ++sh_p1_vrows;
            sh_p1_vcells += vend - vlo;
        }
#endif
        for (col = LB[row]; col <= RB[row]; ++col) {
#if YAMA_SIMD
            if (col == vlo) {
                n = ua[row][0];
                cd_avx2(vlo, vend, pC, pD, pI, cC, cD, fcd, cs, stride, nsym,
                        rt + row, wA + (size_t)row*nsym, (score_t)(n*L*gap_extend));
                if ((col = vend) > RB[row])
                    break;
            }
#endif
            // compute C and flag_c from (row-1,col-1)
            if (col > LB[row-1]) {
                x = pC[col-1];
                y = pD[col-1];
                z = pI[col-1];
                if (col > 1) {	// no gap-open penalty at start
                    if (row > 1 && col > LB[row-2] + 1)
                        x -= DOT4(gb[col], rt[row].Cx);
                    if (col > LB[row-1] + 1)
                        z -= DOT4(gb[col], rt[row].Cz);
                    if (row > 1)
                        y -= vb[col][0]*rt[row].Cy[0] + vb[col][1]*rt[row].Cy[1];
                }
                BEST3(x, y, z, C, flag_c);
                for (k = 0; k < nsym; ++k)	// a zero count adds zero: no branch
                    C += hist[(size_t)col*nsym + k]*wA[(size_t)row*nsym + k];
            } else {
                C = MININT;
                flag_c = 0;
            }

            // compute D and flag_d from (row-1,col)
            x = pC[col];
            y = pD[col];
            z = pI[col];
            if (0 < col && col < N) {	// not for end-gaps
                // here v is 1 and t is B[col]'s own dash flag
                if (row > 1 && col > LB[row-2])
                    x -= vb[col][0]*rt[row].Dx[0] + vb[col][1]*rt[row].Dx[1];
                if (row > 1)
                    y -= rt[row].Dy;
                if (col > LB[row-1])
                    z -= vb[col][0]*rt[row].Dz[0] + vb[col][1]*rt[row].Dz[1];
            }
            BEST3(x, y, z, D, flag_d);
            // Each of n non-dashes in the deleted column
            // is aligned to L dashes.
            n = ua[row][0];	// non-dashes in the deleted column
            D -= n*L*gap_extend;

            cC[col] = C;
            cD[col] = D;
            fcd[col] = flag_c | (flag_d<<2);
        }
#ifdef YAMA_SHADOW
        vlo = vend = -1;
        for (col = LB[row]; col <= RB[row]; ++col) {	// copy of the scalar pass 1
            if (col > LB[row-1]) {
                sx = pC[col-1];
                sy = pD[col-1];
                sz = pI[col-1];
                if (col > 1) {
                    if (row > 1 && col > LB[row-2] + 1)
                        sx -= DOT4(gb[col], rt[row].Cx);
                    if (col > LB[row-1] + 1)
                        sz -= DOT4(gb[col], rt[row].Cz);
                    if (row > 1)
                        sy -= vb[col][0]*rt[row].Cy[0] + vb[col][1]*rt[row].Cy[1];
                }
                BEST3(sx, sy, sz, sV, sf);
                for (sk = 0; sk < nsym; ++sk)
                    sV += hist[(size_t)col*nsym + sk]*wA[(size_t)row*nsym + sk];
            } else {
                sV = MININT;
                sf = 0;
            }
            sC[col] = sV;
            sx = pC[col];
            sy = pD[col];
            sz = pI[col];
            if (0 < col && col < N) {
                if (row > 1 && col > LB[row-2])
                    sx -= vb[col][0]*rt[row].Dx[0] + vb[col][1]*rt[row].Dx[1];
                if (row > 1)
                    sy -= rt[row].Dy;
                if (col > LB[row-1])
                    sz -= vb[col][0]*rt[row].Dz[0] + vb[col][1]*rt[row].Dz[1];
            }
            BEST3(sx, sy, sz, sV, sg);
            n = ua[row][0];
            sV -= n*L*gap_extend;
            sD[col] = sV;
            sfcd[col] = sf | (sg<<2);
            SH_CMP("C", cC[col], sC[col]);
            SH_CMP("D", cD[col], sD[col]);
            SH_CMP("flag_c", fcd[col] & 3, sf);
            SH_CMP("flag_d", (fcd[col] >> 2) & 3, sg);
            SH_CMP("fcd_high_bits", fcd[col] >> 4, 0);
            sh_val(cC[col]);
            sh_val(cD[col]);
        }
#endif

        /* Pass 2: I, which reads C, D and I at (row,col-1), left to right,
        *  then the traceback byte.
        */
        I = MININT;
        flag_i = 0;	// unused
        cI[LB[row]] = I;
        *(tbp++) = fcd[LB[row]] | (flag_i << 4);
        vlo = vend = -1;	// no vector segment
#if YAMA_SIMD && YAMA_SIMD_I
        if (simd && row < M && rowok[row]) {
            vlo = LB[row] + 2;	// both guards below hold from here
            vend = vlo + ((RB[row] - vlo + 1) & ~3);
            if (vend <= vlo)
                vlo = vend = -1;
        }
#endif
#ifdef YAMA_SHADOW
        if (vlo >= 0) {
            ++sh_p2_vrows;
            sh_p2_vcells += vend - vlo;
        }
#endif
        for (col = LB[row] + 1; col <= RB[row]; ++col) {
#if YAMA_SIMD && YAMA_SIMD_I
            if (col == vlo) {
                I = i_avx2(vlo, vend, cC, cD, cI, fcd, tbp, cs, stride, Izc,
                           eK, rt + row, I);
                tbp += vend - vlo;
                if ((col = vend) > RB[row])
                    break;
            }
#endif
            x = cC[col-1];
            y = cD[col-1];
            z = I;
            /* The main complexity is for the quasi-natural
            *  gap costs, which depend on the types (C,D,I)
            *  of the last two edges in a path, i.e., on
            *  the types of the last two nodes.
            */
            if (row < M) {	// not for end-gaps
                if (col > LB[row-1]+1)
                    x -= DOT4(gb[col], rt[row].Ix);
                if (col > LB[row]+1)
                    z -= Izc[col];
                y -= vb[col][0]*rt[row].Iy[0] + vb[col][1]*rt[row].Iy[1];
            }
            /* Precede the insertion edge with a substitution edge
            *  (from the C node at (row-1,col)) if that is best,
            *  else with a deletion edge, else with an insertion.
            */
            BEST3(x, y, z, I, flag_i);
            /* Each of n non-dashes in the inserted column
            * is aligned to K dashes.
            */
            n = vb[col][0];	// non-dashes in the inserted column
            I -= n*K*gap_extend;
            cI[col] = I;

            // pack three flags into one byte of traceback memory
            *(tbp++) = fcd[col] | (flag_i << 4);
        }
#ifdef YAMA_SHADOW
        {   // copy of the scalar pass 2, against the stored I and traceback bytes
            sV = MININT;
            sf = 0;
            col = LB[row];
            SH_CMP("I", cI[col], sV);
            SH_CMP("tback", tback_row[row][col], sfcd[col]);
            sh_val(cI[col]);
            ++sh_cells;
            for (col = LB[row] + 1; col <= RB[row]; ++col) {
                sx = cC[col-1];
                sy = cD[col-1];
                sz = sV;
                if (row < M) {
                    if (col > LB[row-1]+1)
                        sx -= DOT4(gb[col], rt[row].Ix);
                    if (col > LB[row]+1)
                        sz -= Izc[col];
                    sy -= vb[col][0]*rt[row].Iy[0] + vb[col][1]*rt[row].Iy[1];
                }
                BEST3(sx, sy, sz, sV, sf);
                n = vb[col][0];
                sV -= n*K*gap_extend;
                SH_CMP("I", cI[col], sV);
                SH_CMP("tback", tback_row[row][col], sfcd[col] | (sf << 4));
                sh_val(cI[col]);
                ++sh_cells;
            }
            if (tbp != tback_row[row] + RB[row] + 1)
                sh_fail("tbp_advance", row, RB[row], (long long)(tbp - tback_row[row]),
                        RB[row] + 1, LB[row], RB[row], M, N, K, L, sh_p1lo, sh_p1end,
                        vlo, vend);
        }
#endif
        C = cC[RB[row]];
        D = cD[RB[row]];

        sw = pC; pC = cC; cC = sw;
        sw = pD; pD = cD; cD = sw;
        sw = pI; pI = cI; cI = sw;
    }

#ifdef YAMA_TBHASH
    {   /* test instrument: FNV-1a of every traceback byte, and the end scores */
        unsigned long long h_ = 1469598103934665603ULL; size_t q_;
        for (q_ = 0; q_ < tback_size; ++q_) { h_ ^= tback[q_]; h_ *= 1099511628211ULL; }
        printf("TB %016llx %lld %lld %lld\n", h_, C, D, I);
    }
#endif
    script = ckalloc((M+N)*sizeof(uchar));
    nedit = 0;
    // start with best node at grid point (M,N)
    row = M;
    col = N;
    if (C >= D && C >= I)
        node = FLAG_C;
    else if (D >= I)
        node = FLAG_D;
    else
        node = FLAG_I;
    while (row > 0 || col > 0) {
        /* Trace back from a C, D or I node at grid point (row, col).
        *  The type (C, D or I) tells the previous grid point on an
        *  optimal path.  The path's node type at that previous grid
        *  point is packed into the traceback byte.
        */
        if (row < 0 || col < 0) {
            fatal("Error generating edit script.");
        }
        st = tback_row[row][col];	// three flags in one byte
        script[nedit++] = node;
        if (node == FLAG_I) {
            col--;
            node = (st>>4);
        } else if (node == FLAG_D) {
            row--;
            node = (st>>2) & SELECT_CID;
        } else if (node == FLAG_C) {
            row--;
            col--;
            node = st & SELECT_CID;
        } else
            fatal("illegal node type in traceback");
    }

    *OM = m_new = nedit;
    al_new = (uchar **)ckalloc(m_new*sizeof(uchar *)) - 1;
    al_new[1] = (uchar *)ckalloc((size_t)m_new*(K+L)*sizeof(uchar));
    for (i = 2; i <= m_new; ++i)
        al_new[i] = al_new[i-1] + (K+L);

    i = j = m = 0;
    while (--nedit >= 0) {
        if ((k = script[nedit]) == FLAG_C)
            new_col(A[++i], K, B[++j], L, al_new[++m]);
        else if (k == FLAG_I)
            new_col(dashes, K, B[++j], L, al_new[++m]);
        else if (k == FLAG_D)
            new_col(A[++i], K, dashes, L, al_new[++m]);
        else
            fatalf("Illegal edit op: %d", k);
    }
    if (i != M || j != N || m  != m_new)
        fatalf("new_align: i=%d, j=%d, m=%d, M=%d, N=%d, M_new=%d\n",
               i, j, j, M, N, m_new);
    *OAL = al_new;

    free(gb);
    free(vb);
    free(hist);
    free(cnt);
    free(ua);
    free(rt);
    free(wA);
    free(Izc);
    free(tback_row);
    free(tback);
    free(dpbuf);
#if YAMA_SIMD
    free(cs);
    free(rowok);
#endif
#if YAMA_SIMD && YAMA_SIMD_I
    free(eK);
#endif
    free(fcd);
#ifdef YAMA_SHADOW
    free(sC);
    free(sfcd);
#endif
    free(dashes);
    free(script);
}
