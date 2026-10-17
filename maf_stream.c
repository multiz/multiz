/* maf_stream.c -- read a MAF file one reference sequence at a time, without loading it.
 *
 * mafIndexBuild() reads the file once, echoing comment lines exactly as mafReadAll() does in verbose
 * mode, and records each block's byte offset and the src of its top row. mafStream then hands out
 * the blocks of one reference sequence in file order, parsing each with mafNext() only when asked.
 * Memory is about 16 bytes per block plus the blocks the caller currently holds.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "maf.h"
#include "maf_stream.h"

#define HASH_SIZE 65536

static unsigned hash_name(char *s) {
    unsigned h = 5381;
    for (; *s; s++)
        h = h * 33 + (unsigned char)*s;
    return h & (HASH_SIZE - 1);
}

struct mafIndex *mafIndexBuild(char *fileName, int verbose) {
    struct mafIndex *ix = ckalloc(sizeof(struct mafIndex));
    FILE *fp;
    char *line = NULL, src[500];
    size_t cap = 0;
    ssize_t len;
    long here;
    int inblock = 0, wantTop = 0, capBlk = 0, capChr = 0, c, *head, *link, i;

    ix->mf = mafOpen(fileName, 0);      // checks the header line; positioned after it
    fp = ix->mf->fp;
    ix->off = NULL;
    ix->chr = NULL;
    ix->nblocks = ix->nchr = 0;
    ix->chrName = NULL;
    ix->chrN = NULL;
    head = ckalloc(HASH_SIZE * sizeof(int));
    for (i = 0; i < HASH_SIZE; i++)
        head[i] = -1;
    link = NULL;

    for (here = ftell(fp); (len = getline(&line, &cap, fp)) != -1; here += len) {
        if (line[0] == '#') {
            if (verbose && strstr(line, "eof") == NULL)
                printf("%s", line);
            continue;
        }
        if (line[0] == 'a') {
            if (ix->nblocks == capBlk) {
                capBlk = capBlk ? 2 * capBlk : 1024;
                ix->off = realloc(ix->off, capBlk * sizeof(long));
                ix->chr = realloc(ix->chr, capBlk * sizeof(int));
                if (ix->off == NULL || ix->chr == NULL)
                    fatal("out of memory indexing maf");
            }
            ix->off[ix->nblocks] = here;
            ix->chr[ix->nblocks] = -1;
            ix->nblocks++;
            inblock = wantTop = 1;
        } else if (line[0] == '\n' || line[0] == ' ')
            inblock = 0;
        else if (inblock && wantTop && line[0] == 's') {
            if (sscanf(line, "s %499s", src) != 1)
                fatalf("bad component in file %s:\n%s", fileName, line);
            wantTop = 0;
            for (c = head[hash_name(src)]; c >= 0; c = link[c])
                if (strcmp(ix->chrName[c], src) == 0)
                    break;
            if (c < 0) {
                if (ix->nchr == capChr) {
                    capChr = capChr ? 2 * capChr : 64;
                    ix->chrName = realloc(ix->chrName, capChr * sizeof(char *));
                    ix->chrN = realloc(ix->chrN, capChr * sizeof(int));
                    link = realloc(link, capChr * sizeof(int));
                    if (ix->chrName == NULL || ix->chrN == NULL || link == NULL)
                        fatal("out of memory indexing maf");
                }
                c = ix->nchr++;
                ix->chrName[c] = copy_string(src);
                ix->chrN[c] = 0;
                link[c] = head[hash_name(src)];
                head[hash_name(src)] = c;
            }
            ix->chr[ix->nblocks - 1] = c;
            ix->chrN[c]++;
        }
    }
    free(line);
    free(head);
    free(link);

    // a block with no 's' row would make mafNext() fail later; fail now, as mafReadAll() would
    for (i = 0; i < ix->nblocks; i++)
        if (ix->chr[i] < 0)
            fatalf("block with no component in file %s", fileName);

    ix->chrBlk = ckalloc((ix->nchr ? ix->nchr : 1) * sizeof(int *));
    for (c = 0; c < ix->nchr; c++) {
        ix->chrBlk[c] = ckalloc((ix->chrN[c] ? ix->chrN[c] : 1) * sizeof(int));
        ix->chrN[c] = 0;
    }
    for (i = 0; i < ix->nblocks; i++) {
        c = ix->chr[i];
        ix->chrBlk[c][ix->chrN[c]++] = i;
    }
    ix->used = ckalloc(ix->nblocks ? ix->nblocks : 1);
    memset(ix->used, 0, ix->nblocks ? ix->nblocks : 1);
    ix->left = ix->nblocks;
    return ix;
}

int mafIndexChr(struct mafIndex *ix, char *name) {
    int c;

    for (c = 0; c < ix->nchr; c++)
        if (strcmp(ix->chrName[c], name) == 0)
            return c;
    return -1;
}

int mafIndexFirstUnused(struct mafIndex *ix) {
    int i;

    for (i = 0; i < ix->nblocks; i++)
        if (!ix->used[i])
            return i;
    return -1;
}

struct mafAli *mafIndexRead(struct mafIndex *ix, int blk) {
    struct mafAli *a;

    if (fseek(ix->mf->fp, ix->off[blk], SEEK_SET) != 0)
        fatalf("cannot seek in %s", ix->mf->fileName);
    a = mafNext(ix->mf);
    if (a == NULL)
        fatalf("cannot re-read block %d of %s", blk, ix->mf->fileName);
    if (!ix->used[blk]) {
        ix->used[blk] = 1;
        ix->left--;
    }
    return a;
}

void mafStreamOpen(struct mafStream *s, struct mafIndex *ix, int chr) {
    s->ix = ix;
    s->chr = chr;
    s->pos = 0;
    s->lastStart = -1;
}

struct mafAli *mafStreamNext(struct mafStream *s) {
    struct mafAli *a;

    if (s->chr < 0 || s->pos >= s->ix->chrN[s->chr])
        return NULL;
    a = mafIndexRead(s->ix, s->ix->chrBlk[s->chr][s->pos++]);
    if (a->components->start < s->lastStart)
        fprintf(stderr, "warning: %s is not sorted: %s block at %d follows one at %d; "
                "multiz assumes blocks sorted by start and will mis-merge them\n",
                s->ix->mf->fileName, a->components->src, a->components->start, s->lastStart);
    s->lastStart = a->components->start;
    return a;
}

void mafIndexFree(struct mafIndex **pix) {
    struct mafIndex *ix = *pix;
    int c;

    if (ix == NULL)
        return;
    for (c = 0; c < ix->nchr; c++) {
        free(ix->chrName[c]);
        free(ix->chrBlk[c]);
    }
    free(ix->chrName);
    free(ix->chrBlk);
    free(ix->chrN);
    free(ix->off);
    free(ix->chr);
    free(ix->used);
    mafFileFree(&ix->mf);
    free(ix);
    *pix = NULL;
}
