/* yama fuzz driver: reads cases from stdin and prints yama's alignment.
 * A case is "K M", M columns of A (K chars each), "L N", N columns of B,
 * then M+1 values of LB and M+1 values of RB.  Output: OM, then the OM
 * columns, A part then B part.  Built against any multiz tree, so two
 * trees' yama can be compared byte for byte (see run.sh).
 */
#include "util.h"
#include "mz_scores.h"
#include "mz_yama.h"
int main(void) {
    int K, M, L, N, i, OM, *LB, *RB;
    uchar **A, **B, **AL;
    static char buf[1 << 16];
    init_scores70();
    {   /* YAMA_FUZZ_GOP=n replaces every nonzero gap-open entry with n, so with
        *  n near 2^28 the per-row coefficients pass 32 bits and yama must take
        *  its scalar fallback (FITS32) on those rows.  YAMA_FUZZ_GOP_IDX="i,j"
        *  sets only gop[i], gop[j] (zero or not) to n instead.  Under the
        *  shipped table every FITS32 check is implied by another (all nonzero
        *  entries are equal), so a lone check can only be tested with a table
        *  multiz never builds: index 2, GAP(0,0,1,0), enlarges only Ix on rows
        *  of A without dashes; 0 only Cx; 8 only Cx and Cz. */
        const char *g = getenv("YAMA_FUZZ_GOP"), *gi = getenv("YAMA_FUZZ_GOP_IDX");
        if (g && gi) {
            char *e;
            for (const char *q = gi; *q; q = *e ? e + 1 : e) {
                long ix = strtol(q, &e, 10);
                if (e == q || ix < 0 || ix > 15)
                    fatalf("bad YAMA_FUZZ_GOP_IDX %s", gi);
                gop[ix] = atoi(g);
            }
        } else if (g)
            for (i = 0; i < 16; i++)
                if (gop[i])
                    gop[i] = atoi(g);
    }
    {   /* YAMA_FUZZ_SS=f multiplies every substitution score by f.  With f near 2^23 a row's summed scores (wA in
        *  yama) pass 32 bits while the value-bound check still lets the call take the AVX2 path, so the wA FITS32
        *  check is the only thing sending those rows to scalar: the one case that exercises it (11.10 review). */
        const char *f = getenv("YAMA_FUZZ_SS");
        if (f) {
            long long m = atoll(f);
            int a, b;
            for (a = 0; a < SS_NBYTES; a++)
                for (b = 0; b < SS_NBYTES; b++) {
                    long long v = (long long)ss[a][b] * m;
                    if (v > INT_MAX || v < INT_MIN)
                        fatalf("YAMA_FUZZ_SS=%s: a score leaves int", f);
                    ss[a][b] = (int)v;
                }
        }
    }
    while (scanf("%d %d", &K, &M) == 2) {
        A = (uchar **)ckalloc(M*sizeof(uchar *)) - 1;
        for (i = 1; i <= M; i++) { if (scanf("%s", buf) != 1) return 1; A[i] = (uchar *)copy_string(buf); }
        if (scanf("%d %d", &L, &N) != 2) return 1;
        B = (uchar **)ckalloc(N*sizeof(uchar *)) - 1;
        for (i = 1; i <= N; i++) { if (scanf("%s", buf) != 1) return 1; B[i] = (uchar *)copy_string(buf); }
        LB = ckalloc((M+1)*sizeof(int)); RB = ckalloc((M+1)*sizeof(int));
        for (i = 0; i <= M; i++) if (scanf("%d", &LB[i]) != 1) return 1;
        for (i = 0; i <= M; i++) if (scanf("%d", &RB[i]) != 1) return 1;
        yama(A, K, M, B, L, N, LB, RB, &AL, &OM);
        printf("%d\n", OM);
        for (i = 1; i <= OM; i++) printf("%.*s\n", K+L, AL[i]);
        fflush(stdout);
        free(AL[1]); free(AL+1);
        for (i = 1; i <= M; i++) free(A[i]);
        for (i = 1; i <= N; i++) free(B[i]);
        free(A+1); free(B+1); free(LB); free(RB);
    }
    return 0;
}
