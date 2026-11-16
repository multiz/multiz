#ifndef MZ_SCORES_H
#define MZ_SCORES_H

#include "maf.h"

// defined once, in mz_scores.c (a definition here is one per includer, which
// gcc 10 and later reject at link time)
extern int **ss,  // SS_NBYTES x SS_NBYTES array of substitution scores
       *gop;      // 16-position array of quasi-natural gap-open penalties

extern int gap_open, gap_extend;

/* ss covers every byte value: the MAF parser accepts any non-blank byte, and a byte of 0x80 or above used to
*  index past a 128-row table and crash multiz (11.10 code review, finding 14). Rows and columns 128-255 hold the
*  same values as any other non-ACGT character, so such a byte scores as N does and ASCII input scores as before.
*  SS_NCHARS stays 128: yama's symbol-class compare and the |ss| maximum loop over it, and both are exact over
*  128, because row and column 128-255 entries repeat values the first 128 already hold (filler, or the gap
*  extension against '-'), so two columns that agree on rows 0-127 agree on all 256. */
#define SS_NBYTES 256	// ss is allocated and filled SS_NBYTES x SS_NBYTES
#define SS_NCHARS 128	// the rows yama compares and scans; see above
#define SS(c,d) ss[c][d]
#define GAP(s,t,u,v) gop[(s<<3)+(t<<2)+(u<<1)+v]
#define GAP2(s,t,u,v) GAP((s == '-'), (t == '-'), (u == '-'), (v == '-'))

void init_scores70();
void init_scores85();
double mafScoreRange(struct mafAli *maf, int start, int size);

#endif
