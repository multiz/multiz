#define _POSIX_C_SOURCE 200809L  /* popen, fileno, pid_t: hidden by strict ISO C */
/*
 *  blastzWrapper.c version 10
 *
 *  -- wrapping of blastz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "multi_util.h"
#include "seq.h"

#define BZ "lastz"    // we use lastz instead of blastz
#define VERSION 16

// replace the first argument with the new first_fn; newcmd holds n bytes, and a longer command is an error
static char* modify_cmdline(int argc, char** argv, char* first_fn, char* second_fn, char* newcmd, size_t n) {
    size_t len;
    int i;
    cksnprintf(newcmd, n, "%s %s %s", BZ, first_fn, second_fn);
    for (i=3; i<argc; i++) {
        len = strlen(newcmd);
        cksnprintf(newcmd + len, n - len, " %s", argv[i]);
    }
    return newcmd;
}

/* lastz's exit status on the per-contig path, which pipes it through grep: the shell appends a "#:status N"
*  line (grep keeps it, as it holds no "eof"), so pclose's status, grep's, no longer hides lastz's failure.
*  -1: no status line arrived. */
static int bz_status = -1;

// switch the first and second sequence
void replace_reverse_bz(FILE* fpr, FILE* fpw, char* replace_str, int contig, int rev) {
    char buf[1000], buf1[1000], buf2[1000], cmd[200], spe1[200], spe2[200];
    int pct, dir;
    maf_pos_t b1, b2, e1, e2, beg, end;

    while ( fgets(buf, 1000, fpr) != NULL) {
        if ( rev==1 && strncmp(buf, "d {", 3) == 0) {
            fprintf(fpw, "%s", buf);
            ckfgets(buf, 1000, fpr, "lav input");
            sscanf(buf, "  \"%s %s %s", cmd, spe1, spe2);
            fprintf(fpw, "  \"%s %s %s\n", cmd, spe2, spe1);
        } else if ( rev==1 && strncmp(buf, "h {",3)==0) {
            fprintf(fpw, "%s", buf);
            ckfgets(buf1, 1000, fpr, "lav input");
            ckfgets(buf2, 1000, fpr, "lav input");
            fprintf(fpw, "%s", buf2);
            fprintf(fpw, "%s", buf1);
        } else if ( strncmp(buf, "s {", 3) == 0) {
            fprintf(fpw, "%s", buf);
            ckfgets(buf1, 1000, fpr, "lav input");
            sscanf(buf1, "%*s %lld %lld %d %*d", &beg, &end, &dir);
            sprintf(buf1, "  \"%s\" %lld %lld %d %d\n", replace_str, beg, end, dir, contig);
            ckfgets(buf2, 1000, fpr, "lav input");
            if (rev==1) {
                fprintf(fpw, "%s", buf2);
                fprintf(fpw, "%s", buf1);
            } else {
                fprintf(fpw, "%s", buf1);
                fprintf(fpw, "%s", buf2);
            }
        } else if ( rev==1 && strncmp(buf, "a {", 3)==0) {
            fprintf(fpw, "%s", buf);
            ckfgets(buf, 1000, fpr, "lav input");
            fprintf(fpw, "%s", buf);  // s
            ckfgets(buf, 1000, fpr, "lav input");    // b
            sscanf(buf, "  b %lld %lld", &b1, &b2);
            fprintf(fpw, "  b %lld %lld\n", b2, b1);
            ckfgets(buf, 1000, fpr, "lav input");    // e
            sscanf(buf, "  e %lld %lld", &e1, &e2);
            fprintf(fpw, "  e %lld %lld\n", e2, e1);
            ckfgets(buf, 1000, fpr, "lav input");
            while ( strncmp(buf, "  l", 3)==0) {
                sscanf(buf, "  l %lld %lld %lld %lld %d", &b1, &b2, &e1, &e2, &pct);
                fprintf(fpw, "  l %lld %lld %lld %lld %d\n", b2, b1, e2, e1, pct);
                ckfgets(buf, 1000, fpr, "lav input");
            }
            fprintf(fpw, "%s", buf);
        } else if ( strncmp(buf, "#:status ", 9) == 0)
            bz_status = atoi(buf + 9);
        else if ( strncmp(buf, "#:eof", 5) != 0)
            fprintf(fpw, "%s", buf);
    }
}

int main(int argc, char** argv) {
    char name[64], cmdline[8192], inner[8192], seq_file1[500], seq_file2[500];
    int c1, c2, reverse=0, tmp, fd, contig=1;
    SEQ* sf;
    FILE* pf, *fpw;

    sprintf(name, "blastzWrapper.v%d", VERSION);
    argv0 = name;   // its own buffer: argv0 was cmdline, which then held each lastz command

    if ( argc < 3 )
        fatal(" -- wrapper of blastz, passing all arguments to blastz.\nargs: seqfile1 seqfile2 [options]");

    sf = seq_open(argv[1]);
    for (c1=0; seq_read(sf)!=0; c1++)
        ;
    seq_close(sf);
    sf = seq_open(argv[2]);
    for (c2=0; seq_read(sf)!=0; c2++)
        ;
    seq_close(sf);

    if ( c1 > c2) {
        reverse = 1;
        cksnprintf(seq_file1, sizeof seq_file1, "%s", argv[2]);
        cksnprintf(seq_file2, sizeof seq_file2, "%s", argv[1]);
        tmp = c1;
        c1 = c2;
        c2 = tmp;
    } else {
        cksnprintf(seq_file1, sizeof seq_file1, "%s", argv[1]);
        cksnprintf(seq_file2, sizeof seq_file2, "%s", argv[2]);
    }

    if (c1 == 1) {
        modify_cmdline(argc, argv, seq_file1, seq_file2, cmdline, sizeof cmdline);
        if (reverse == 0)
            return system(cmdline) == 0 ? 0 : 1;   // lastz's failure is blastzWrapper's
        if ((pf = popen(cmdline, "r")) == NULL)
            fatalf("cannot run '%s'", cmdline);
        replace_reverse_bz(pf, stdout, seq_file1, contig, 1);
        if (pclose(pf) != 0)    // lastz's own status; ignored, it ended the output with "#:eof" as if it had worked
            fatalf("command '%s' failed", cmdline);
    } else {
        sf = seq_open(seq_file1);
        while ( seq_read(sf)) {
            fpw = tmpfile();
            fprintf(fpw, "%s\n", sf->header);
            fprintf(fpw, "%s\n", sf->seq);
            fseek(fpw, 0, SEEK_SET);
            fd = fileno(fpw);
            dup2(fd, 0);
            modify_cmdline(argc, argv, "/dev/stdin", seq_file2, inner, sizeof inner);
            cksnprintf(cmdline, sizeof cmdline, "{ %s; echo \"#:status $?\"; } | grep -v eof", inner);
            if ((pf = popen(cmdline, "r")) == NULL)
                fatalf("cannot run '%s'", cmdline);
            bz_status = -1;
            if ( reverse == 1)
                replace_reverse_bz(pf, stdout, seq_file1, contig, 1);
            else
                replace_reverse_bz(pf, stdout, seq_file1, contig, 0);
            if (pclose(pf) != 0 || bz_status != 0)
                fatalf("command '%s' failed on sequence %d of %s", inner, contig, seq_file1);
            fclose(fpw);
            contig++;
        }
        seq_close(sf);
    }

    fprintf(stdout, "#:eof\n");
    return 0;
}
