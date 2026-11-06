#ifndef MAF_H
#define MAF_H
/*
 *   maf.h version 12
 */
// mini-maf.h - Multiple alignment format, header file stolen from Jim Kent and
// abused
#include <stdio.h>

#include <limits.h>
#define MAX_INT INT_MAX
#define MIN_INT INT_MIN

/* A sequence coordinate or length. 64-bit so a chromosome may exceed 2^31 - 1 bases (11.6). */
typedef long long maf_pos_t;

/* Three-way comparison of two coordinates, for qsort comparators: a difference would not fit in an int. */
static inline int maf_pos_cmp(maf_pos_t a, maf_pos_t b) {
    return (a > b) - (a < b);
}

struct mafFile
            /* A file full of multiple alignments. */
{
    struct mafFile *next;
    int version;		// Required
    char *scoring;		// Name of scoring scheme.
    struct mafAli *alignments;	// Possibly empty list of alignments.
    char *fileName;
    int line_nbr;
    int verbose;
    FILE *fp; 			// Open file if any. NULL except while parsing.
    char *pending;		// an 'a' line read while finishing the previous block
};

void mafFileFree(struct mafFile **pObj);
/* Free up a maf file including closing file handle if necessary. */

struct mafAli
            /* A multiple alignment. */
{
    struct mafAli *next;
    double score;
    struct mafComp *components;	/* List of components of alignment */
    int textSize;	 /* Size of text in each component. */
    int chain_len;
};

void mafAliFree(struct mafAli **pObj);
/* Free up a maf alignment. */

struct mafComp
            /* A component of a multiple alignment. */
{
    struct mafComp *next;
    char *name;        /* comman name of sequence source. */
    char *src;	 /* Name of sequence source.  */
    char *text;        /* The sequence including dashes. */
    char* contig;
    maf_pos_t* mafPosMap;
    maf_pos_t srcSize;       /* Size of sequence source.  */
    maf_pos_t start;	 /* Start within sequence. Zero based. If strand is - is relative to src end. */
    maf_pos_t size;	         /* Size in sequence (does not include dashes).  */
    short nameID;
    char strand;       /* Strand of sequence.  Either + or -*/
    char paralog;
};

void mafCompFree(struct mafComp **pObj);
/* Free up a maf component. */

struct mafFile *mafOpen(char *fileName, int verbose);
/* Open up a .maf file for reading.  Read header and
 * verify. Prepare for subsequent calls to mafNext().
 * Prints error message and aborts if there's a problem. */

struct mafAli *mafNext(struct mafFile *mafFile);
/* Return next alignment in file or NULL if at end.
 * This will close the open file handle at end as well. */

struct mafFile *mafReadAll(char *fileName, int verbose);
/* Read in full maf file */

void mafWriteStart(FILE *f, char *scoring);
/* Write maf header and scoring scheme name (may be null) */

void mafWriteEnd(FILE *f);

void mafWrite(FILE *f, struct mafAli *maf);
/* Write next alignment to file. */

struct mafAli* mafColDashRm(struct mafAli *a);
/* remove columns which contain all dashes. */

struct mafAli* mafRowDashRm(struct mafAli *a);
/* remove components which contain all dashes. */

struct mafComp *mafNewComp(char *src, maf_pos_t start, maf_pos_t size, char strand, maf_pos_t srcSize, int len, char paralog, char* name, char* contig);
/* construct a new component. */

struct mafAli* make_part_ali(struct mafAli* template, int cbeg, int cend);

struct mafComp *mafCpyComp(struct mafComp*);
struct mafAli *mafNewAli(double score, int textSize);
struct mafAli* duplicate_ali(struct mafAli* template);

#endif
