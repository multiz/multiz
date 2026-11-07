#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include "speciesTree.h"
#include "util.h"
#include "multi_util.h"

extern int force;
extern int execute;
extern int verbose;
extern char **environ;

extern char* PREFIX;
extern char* OPERAT;

#define BIG_BUF_SIZE 50000
#define NAMESIZE_ST 500	// matches roast's NAMESIZE: temp-file names are built in buffers this size

// do_cmd -- print and possible execute a command
void do_cmd(const char *fmt, ...) {
    char *buf;
    va_list ap;
    int n;

    // size the command first: a fixed buffer overflowed silently on long
    // command lines (many input files, or the echoed header of a long run)
    va_start(ap, fmt);
    n = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (n < 0)
        fatal("do_cmd: cannot format command");
    buf = ckalloc((size_t)n + 1);
    va_start(ap, fmt);
    (void)vsnprintf(buf, (size_t)n + 1, fmt, ap);
    va_end(ap);

    if (verbose) {
        (void)printf("%s\n", buf);
        (void)fflush(stdout);
    }
    if (execute)
        if (system(buf) != 0)
            if (force == 0)
                fatalf("command '%s' failed", buf);
    free(buf);
}


/* ------------- native versions of the shell file steps (see speciesTree.h) ------------- */

static void native_note(const char *cmd, int ok);
static void native_show(const char *cmd);

void native_run(const char *words, const char *const *args, const char *out, int append) {
    char *cmd, *copy, *tok, *save = NULL, **argv;
    size_t len, nw = 0, na = 0, i;
    posix_spawn_file_actions_t fa;
    pid_t pid;
    int status = 0, ok;

    for (na = 0; args[na] != NULL; na++)
        ;
    len = strlen(words) + strlen(out) + 8;
    for (i = 0; i < na; i++)
        len += strlen(args[i]) + 1;
    cmd = ckalloc(len);
    strcpy(cmd, words);              // the text do_cmd would have printed and run
    for (i = 0; i < na; i++) {       // one space before each argument, as the "%s %s" formats gave
        strcat(cmd, " ");
        strcat(cmd, args[i]);
    }
    strcat(cmd, append ? " >> " : " > ");
    strcat(cmd, out);
    native_show(cmd);
    if (!execute) {
        free(cmd);
        return;
    }
    copy = copy_string(words);
    argv = ckalloc((strlen(words) + na + 2) * sizeof(char *));
    for (tok = strtok_r(copy, " \t", &save); tok != NULL; tok = strtok_r(NULL, " \t", &save))
        argv[nw++] = tok;
    if (nw == 0)
        fatalf("command '%s' failed", cmd);
    for (i = 0; i < na; i++)
        argv[nw + i] = (char *)args[i];
    argv[nw + na] = NULL;
    fflush(stdout);                  // the child inherits stdout's descriptor, not its buffer
    fflush(stderr);
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_addopen(&fa, 1, out, O_WRONLY | O_CREAT | (append ? O_APPEND : O_TRUNC), 0666);
    ok = posix_spawnp(&pid, argv[0], &fa, NULL, argv, environ) == 0;
    if (!ok)
        fprintf(stderr, "%s: %s\n", argv[0], strerror(errno));
    posix_spawn_file_actions_destroy(&fa);
    if (ok) {
        while (waitpid(pid, &status, 0) == -1)
            if (errno != EINTR) {
                ok = 0;
                break;
            }
        ok = ok && WIFEXITED(status) && WEXITSTATUS(status) == 0;
    }
    native_note(cmd, ok);
    free(argv);
    free(copy);
    free(cmd);
}

static void native_note(const char *cmd, int ok) {
    if (!ok && force == 0)
        fatalf("command '%s' failed", cmd);
}

static void native_show(const char *cmd) {
    if (verbose) {
        (void)printf("%s\n", cmd);
        (void)fflush(stdout);
    }
}

// copy every byte of from to the end of an open stream; returns 0 if from cannot be read
static int append_file(const char *from, FILE *to) {
    FILE *in = fopen(from, "r");
    char *buf;
    size_t n;
    int ok = 1;

    if (in == NULL)
        return 0;
    buf = ckalloc(1 << 20);
    while ((n = fread(buf, 1, 1 << 20, in)) > 0)
        if (fwrite(buf, 1, n, to) != n) {
            ok = 0;
            break;
        }
    if (ferror(in))
        ok = 0;
    free(buf);
    fclose(in);
    return ok;
}

// cp from to: the source is opened first, so a missing source creates no target, as with cp
static int copy_file(const char *from, const char *to) {
    FILE *in, *out;
    int ok;

    if ((in = fopen(from, "r")) == NULL)
        return 0;
    fclose(in);
    if ((out = fopen(to, "w")) == NULL)
        return 0;
    ok = append_file(from, out);
    if (fclose(out) != 0)
        ok = 0;
    return ok;
}

void native_mv(const char *from, const char *to) {
    char cmd[2 * NAMESIZE_ST + 8];
    int ok;

    snprintf(cmd, sizeof cmd, "mv %s %s", from, to);
    native_show(cmd);
    if (!execute)
        return;
    ok = (rename(from, to) == 0);
    if (!ok && errno == EXDEV)      // another file system: copy, then remove, as mv does
        ok = copy_file(from, to) && unlink(from) == 0;
    native_note(cmd, ok);
}

void native_rm2(const char *a, const char *b) {
    char cmd[2 * NAMESIZE_ST + 8];

    snprintf(cmd, sizeof cmd, "rm -f %s %s", a, b);
    native_show(cmd);
    if (!execute)
        return;
    (void)unlink(a);    // rm -f: a missing file is not an error
    (void)unlink(b);
}

void native_cp(const char *from, const char *to) {
    char cmd[2 * NAMESIZE_ST + 8];

    snprintf(cmd, sizeof cmd, "cp %s %s", from, to);
    native_show(cmd);
    if (!execute)
        return;
    native_note(cmd, copy_file(from, to));
}

/* grep -v: copy each line that does not contain pat anywhere, adding a newline to a last
*  line that lacks one. grep exits 1 when it selects no line and 2 when a file cannot be
*  read, after printing "grep: FILE: No such file or directory" (or the system's reason).
*/
void native_grep_v_append(const char *pat, const char *in1, const char *in2, const char *out) {
    char cmd[4 * NAMESIZE_ST + 32];
    const char *ins[2];
    FILE *o, *in;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int k, selected = 0, err = 0;
    char *obuf, *ibuf;

    if (in2 == NULL)
        snprintf(cmd, sizeof cmd, "grep -v %s %s >> %s", pat, in1, out);
    else
        snprintf(cmd, sizeof cmd, "grep -v -h %s %s %s >> %s", pat, in1, in2, out);
    native_show(cmd);
    if (!execute)
        return;
    if ((o = fopen(out, "a")) == NULL) {    // the shell's >> failed: grep never ran
        native_note(cmd, 0);
        return;
    }
    obuf = ckalloc(1 << 20);
    setvbuf(o, obuf, _IOFBF, 1 << 20);
    ibuf = ckalloc(1 << 20);
    ins[0] = in1;
    ins[1] = in2;
    for (k = 0; k < 2 && ins[k] != NULL; k++) {
        if ((in = fopen(ins[k], "r")) == NULL) {
            fprintf(stderr, "grep: %s: %s\n", ins[k], strerror(errno));
            err = 1;
            continue;
        }
        setvbuf(in, ibuf, _IOFBF, 1 << 20);
        while ((len = getline(&line, &cap, in)) != -1) {
            if (strstr(line, pat) != NULL)
                continue;
            selected = 1;
            fwrite(line, 1, (size_t)len, o);
            if (line[len-1] != '\n')
                fputc('\n', o);
        }
        fclose(in);
    }
    free(line);
    if (fclose(o) != 0)
        err = 1;
    free(obuf);
    free(ibuf);
    native_note(cmd, !err && selected);
}

int  parseSpeciesTree(char* treeStr, TreeNode* tree, int nbz, char** bz_file, int* id, char* buf, int (*operation)(TreeNodePtr, TreeNodePtr, int, int, char **)) {
    char name_list[10000], *q, *p;
    char fn_left[NAMESIZE_ST], fn_right[NAMESIZE_ST], fn_out[NAMESIZE_ST], fn_head[NAMESIZE_ST], fn_child[NAMESIZE_ST];
    NameListPtr n;
    int top=-1, i, j;

    for (q = treeStr; *q != '\0'; ++q) {
        if (*q == '(') {
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            tree[top].type = '(';
        } else if (*q == ')') {
            if (top < 1 || tree[top].type != 0 || tree[top-1].type != '(') {
                q[1] = '\0';
                fatalf("parse error: %s", treeStr);
            }
            tree[top-1] = tree[top];
            --top;
        } else if (isalpha(*q)) {
            /* leaf species */
            if (++top >= 1000)
                fatal("parse_tree: stack overflow");
            p = buf;
            while (isalpha(*q) || isdigit(*q) || *q == '_' || *q == '.')
                *p++ = *q++;
            --q;
            *p = '\0';
            n = ckalloc(sizeof(NameList));
            n->name = copy_string(buf);
            n->next = NULL;
            tree[top].id = -1;
            tree[top].type = 0;
            tree[top].names = n;
        } else if (*q != ' ')
            fatalf("improper character in tree specification: %c", *q);
        if (top > 0 && tree[top-1].type == 0 && tree[top].type == 0) {
            // Process an internal node by merging its children.
            // One optimization is to identify nodes where both
            // children are leaves, to avoid a use of multiz.
            i = tree[top-1].id;
            j = tree[top].id;
            snprintf(fn_left, sizeof fn_left, "%sleft.maf%d", PREFIX, *id);
            snprintf(fn_right, sizeof fn_right, "%sright.maf%d", PREFIX, *id);
            snprintf(fn_out, sizeof fn_out, "%s%s%d", PREFIX, OPERAT, *id);
            snprintf(fn_head, sizeof fn_head, "%shead", PREFIX);
            if (i >= 0) {
                snprintf(fn_child, sizeof fn_child, "%s%s%d", PREFIX, OPERAT, i);
                native_mv(fn_child, fn_left);
            }
            if (j >= 0) {
                snprintf(fn_child, sizeof fn_child, "%s%s%d", PREFIX, OPERAT, j);
                native_mv(fn_child, fn_right);
            }
            native_cp(fn_head, fn_out);
            operation(&(tree[top-1]), &(tree[top]), *id, nbz, bz_file);
            if (i >= 0 || j >= 0) {
                if ( execute == 1) {
                    force = 1;
                    native_grep_v_append("maf", fn_left, NULL, fn_out);
                    native_grep_v_append("maf", fn_right, NULL, fn_out);
                    force = 0;
                }
                name_list[0] = '\0';
                for (n = tree[top-1].names; n != NULL; n = n->next) {
                    strcat(name_list, " ");
                    strcat(name_list, n->name);
                }
                strcat(name_list, " +");
                for (n = tree[top].names; n != NULL; n = n->next) {
                    strcat(name_list, " ");
                    strcat(name_list, n->name);
                }
            }
            /* merge the lists of names */
            n = tree[top-1].names;
            if (n == NULL)
                fatal("empty list of names");
            while (n->next != NULL)
                n = n->next;
            n->next = tree[top].names;
            --top;
            tree[top].id = (*id)++;
        }
    }
    return top;
}

/*
NameListPtr formNameList(char* str) {
  char *ptr, tail='a';
  NameListPtr nname, ret=NULL;

  while ( *str == ' ')
    str++;

  while ( (ptr=strchr(str, ' ')) != NULL || (ptr=strchr(str, '\0')) != NULL ) {
    if ( *ptr == '\0' )
      tail = '\0';
    *ptr = '\0';
    nname = (NameListPtr)malloc(sizeof(NameList));
    nname->name = copy_string(str);
    nname->next = ret;
    ret = nname;
    if ( tail == '\0' )
      break;
    *ptr = ' ';
    str = ptr+1;
    while (*str == ' ')
      str++;
  }
  return ret;
}
*/
