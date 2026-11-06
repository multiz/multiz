// mz_scores.c -- handles some details of alignment scores for multiz

#include "util.h"
#include "mz_scores.h"

int **ss, *gop, gap_open, gap_extend;
int **ss70, **ss85, *gop70, *gop85;

static const uchar nchars[] = "ACGT";
static const int HOX70[4][4] = {	// for human-rodent
                                   {  91, -114,  -31, -123 },
                                   {-114,  100, -125,  -31 },
                                   { -31, -125,  100, -114 },
                                   {-123,  -31, -114,   91 },
                               };

static const int HOX85[4][4] = {	// for mouse-rat alignments
                                   {  86, -135,  -68, -157 },
                                   {-135,  100, -148,  -68 },
                                   { -68, -148,  100, -135 },
                                   {-157,  -68, -135,   86 },
                               };

#define GAP_OPEN70	400
#define GAP_EXTEND70	30

#define GAP_OPEN85	600
#define GAP_EXTEND85	50

#define UNSPECIFIED	-100	// score for an unspecified aligned pair

#define CLEN(s) (sizeof((s))-1)
#define NACHARS 128

static void init_scores(const int s[4][4], int **ss, int *gop,  int filler,
                        int gap_op, int gap_ex) {
    int i, j, a, b, A, B, X, D;

    // first, fill in ss, the substitution-score matrix
    for (i = 0; i < NACHARS; ++i)
        for (j = 0; j < NACHARS; ++j)
            ss[i][j] = filler;
    for (i = 0; i < (signed)CLEN(nchars); ++i) {
        A = nchars[i];
        a = tolower(A);
        for (j = 0; j < (signed)CLEN(nchars); ++j) {
            B = nchars[j];
            b = tolower(B);
            ss[A][B] = ss[a][B] = ss[A][b] = ss[a][b] =
                                                 s[i][j];
        }
    }
    for (i = 0; i < NACHARS; ++i)
        ss['-'][i] = ss[i]['-'] = -gap_ex;
    ss['-']['-'] = 0;

    // initialize the "gap-open penalty" array, for quasi-natural gap costs
    for (i = 0; i < 16; ++i)
        gop[i] = 0;
    D = 1;	// dash
    X = 0;	// anything other than '-'
    /* The six gap-open configurations are:
    	xx	x-	x-	-x	-x	--
    	x-	xx	-x	x-	--	-x
    	 1	 2	 3	 4	 5	 6
       The last two may result in a gap being counted twice, as with:
    	xx-xx
    	x---x
    	 1 5
       (A digit indicates the second column of a column-pair that conforms
       to the indicated rule, 1-6.)
       However, it permits a gap to be detected by inspection of adjacent
       columns in a case like:
    	x-----xx
    	x------x
    	      5
       GAP's arguments give column 1, followed by column 2.
    */
    GAP(X,X,X,D) = GAP(X,X,D,X) = GAP(X,D,D,X) =  GAP(D,X,X,D) =
                                      GAP(D,D,X,D) = GAP(D,D,D,X) = gap_op;
    gap_extend = gap_ex;
}

static int **alloc_scores(int **gop_ptr) {
    int **matrix, i;

    *gop_ptr = ckalloc(16*sizeof(int));
    matrix = ckalloc(NACHARS*sizeof(int *));
    matrix[0] = ckalloc(NACHARS*NACHARS*sizeof(int));
    for (i = 1; i < NACHARS; ++i)
        matrix[i] = matrix[i-1] + NACHARS;
    return matrix;
}

void init_scores70() {
    static int init = 1;

    if (init) {
        init = 0;
        ss70 = alloc_scores(&gop70);
        init_scores(HOX70, ss70, gop70, UNSPECIFIED, GAP_OPEN70,
                    GAP_EXTEND70);
    }
    ss = ss70;
    gop = gop70;
    gap_open = GAP_OPEN70;
    gap_extend = GAP_EXTEND70;
}

void init_scores85() {
    static int init = 1;

    if (init) {
        init = 0;
        ss85 = alloc_scores(&gop85);
        init_scores(HOX85, ss85, gop85, UNSPECIFIED, GAP_OPEN85,
                    GAP_EXTEND85);
    }
    ss = ss85;
    gop = gop85;
    gap_open = GAP_OPEN85;
    gap_extend = GAP_EXTEND85;
}

static double mafScoreRange_pairs(struct mafAli *maf, int start, int size) {
    uchar ai, ar, bi, br;
    int i;
    double score;
    struct mafComp *c1, *c2;

    score = 0.0;
    for (i = start; i < start+size; ++i) {
        for (c1 = maf->components; c1 != NULL; c1 = c1->next) {
            br = c1->text[i];
            for (c2 = c1->next; c2 != NULL; c2 = c2->next) {
                bi = c2->text[i];
                score += SS(br,bi);
                if (i > 0) {
                    ar = c1->text[i-1];
                    ai = c2->text[i-1];
                    score -= GAP2(ar,ai,br,bi);
                }
            }
        }
    }
    return score;
}

/* Counting version: per column, count rows per symbol and per (prev dash, cur dash) class. */
static double mafScoreRange_counts(struct mafAli *maf, int start, int size) {
    long long cnt[256], g[4], total = 0;
    int seen[256], nseen, i, a, b, p, q;
    struct mafComp *c;

    for (a = 0; a < 256; ++a) cnt[a] = 0;
    for (i = start; i < start+size; ++i) {
        nseen = 0;
        g[0] = g[1] = g[2] = g[3] = 0;
        for (c = maf->components; c != NULL; c = c->next) {
            uchar br = c->text[i];
            if (cnt[br]++ == 0) seen[nseen++] = br;
            if (i > 0)
                g[((c->text[i-1] == '-') << 1) | (br == '-')]++;
        }
        for (a = 0; a < nseen; ++a) {
            int sa = seen[a];
            total += cnt[sa]*(cnt[sa]-1)/2 * (long long)SS(sa,sa);
            for (b = a+1; b < nseen; ++b)
                total += cnt[sa]*cnt[seen[b]] * (long long)SS(sa,seen[b]);
        }
        if (i > 0)
            for (p = 0; p < 4; ++p) {
                if (g[p] == 0) continue;
                /* class p = (prev dash, cur dash) for row 1, q for row 2; GAP(s,t,u,v): s,t prev, u,v cur */
                total -= g[p]*(g[p]-1)/2 * (long long)GAP((p>>1), (p>>1), (p&1), (p&1));
                for (q = p+1; q < 4; ++q)
                    total -= g[p]*g[q] * (long long)GAP((p>>1), (q>>1), (p&1), (q&1));
            }
        for (a = 0; a < nseen; ++a) cnt[seen[a]] = 0;
    }
    return (double)total;
}

/* The counting version sums unordered row pairs, which equals the pair loop only when the substitution and
   gap tables are symmetric under swapping the two rows. Checked once per table; an asymmetric table falls back to
   the pair loop, so the result never changes. */
static int pairs_needed(void) {
    static int **checked = NULL;
    static int need = 0;
    int a, b, s, t, u, v;
    if (checked == ss)
        return need;
    need = 0;
    for (a = 0; a < NACHARS && !need; ++a)
        for (b = 0; b < a; ++b)
            if (ss[a][b] != ss[b][a]) { need = 1; break; }
    for (s = 0; s < 2 && !need; ++s)
        for (t = 0; t < 2; ++t)
            for (u = 0; u < 2; ++u)
                for (v = 0; v < 2; ++v)
                    if (GAP(s,t,u,v) != GAP(t,s,v,u)) need = 1;
    checked = ss;
    return need;
}

double mafScoreRange(struct mafAli *maf, int start, int size) {
    double s;
    if (start < 0 || size <= 0 || start+size > maf->textSize)
        fatalf("mafScoreRange: start = %d, size = %d, textSize = %d\n",
               start, size, maf->textSize);
    if (ss == NULL)
        fatal("mafScoreRange: scores not initialized");
#ifdef MSR_CHECK
    {
        double s0 = mafScoreRange_pairs(maf, start, size);
        s = mafScoreRange_counts(maf, start, size);
        if (s != s0) fatalf("MSR MISMATCH %.1f %.1f start=%d size=%d", s0, s, start, size);
    }
#else
    s = pairs_needed() ? mafScoreRange_pairs(maf, start, size) : mafScoreRange_counts(maf, start, size);
#endif
    return s;
}
