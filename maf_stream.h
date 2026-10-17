#ifndef MAF_STREAM_H
#define MAF_STREAM_H

#include <stdio.h>
#include "maf.h"

/* An index of a MAF file: the byte offset of every block, grouped by the src of its top row. Lets a
   program visit the blocks of one reference sequence at a time, in file order, without holding the
   file in memory. */
struct mafIndex {
    struct mafFile *mf;     // the file, re-read block by block (verbose off: comments were echoed once)
    long *off;              // offset of each block's 'a' line, in file order
    int *chr;               // reference-sequence id of each block, in file order
    char *used;             // block already handed out by a reader
    int nblocks, left;      // all blocks; blocks not yet handed out
    char **chrName;         // reference-sequence names, by id, in order of first appearance
    int **chrBlk, *chrN;    // per id: block numbers in file order, and their count
    int nchr;
};

/* A reader over the blocks of one reference sequence, in file order. */
struct mafStream {
    struct mafIndex *ix;
    int chr, pos;
    int lastStart;          // for the sort-order check
};

struct mafIndex *mafIndexBuild(char *fileName, int verbose);
int mafIndexChr(struct mafIndex *ix, char *name);   // id of a reference sequence, or -1
int mafIndexFirstUnused(struct mafIndex *ix);       // block number of the first unused block, or -1
struct mafAli *mafIndexRead(struct mafIndex *ix, int blk);
void mafStreamOpen(struct mafStream *s, struct mafIndex *ix, int chr);
struct mafAli *mafStreamNext(struct mafStream *s);
void mafIndexFree(struct mafIndex **pix);

#endif
