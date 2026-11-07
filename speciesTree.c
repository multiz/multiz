#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
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
static int native_shell(const char *cmd);

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
    if (execute) {
        int ok = native_shell(buf);   // -1: not a form it runs; the shell does, as before
        if (ok < 0)
            ok = (system(buf) == 0);
        if (!ok && force == 0)
            fatalf("command '%s' failed", buf);
    }
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

/* The shell step was grep -v, meant to drop the MAF comment lines that name "maf" or
*  "eof" (the ##maf header, the ##eof terminator, the programs' "# maf_project ..."
*  lines). grep -v also dropped any data line containing the pattern, so a sequence
*  named, for example, "xmaf1" lost its rows and left blocks short of a species. Only
*  comment lines (those starting with '#') are dropped now. Otherwise this does what
*  grep did: a newline is added to a last line that lacks one, the step counts as failed
*  when no line is copied or a file cannot be read, and an unreadable file is reported
*  as "grep: FILE: reason".
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
            if (line[0] == '#' && strstr(line, pat) != NULL)
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

/* ---------- native_shell: run the simple commands do_cmd builds without /bin/sh ----------
*  Recognised forms (words separated by spaces or tabs, no quoting or expansion):
*    PROGRAM ARGS... [> FILE | >> FILE]      started with posix_spawnp
*    mv A B          rm -f FILES...          touch FILES...
*    cat FILES... > FILE                     grep -v [-h] eof|maf FILES... > FILE | >> FILE
*  Anything else (quotes, globs, $, pipes, other options) returns -1 and do_cmd uses
*  the shell, as before. Returns 1 when the command succeeded and 0 when the shell
*  version would have exited non-zero; error messages follow the tools' own wording.
*/
#define NS_MAXW 256

static int ns_spawn(char **argv, const char *out, int append) {
    posix_spawn_file_actions_t fa;
    pid_t pid;
    int status = 0, ok;

    fflush(stdout);
    fflush(stderr);
    posix_spawn_file_actions_init(&fa);
    if (out != NULL)
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
    return ok;
}

// grep -v PAT FILES...: comment lines containing PAT go (see native_grep_v_append)
static int ns_grep_v(const char *pat, char **in, int nin, const char *out, int append) {
    FILE *o, *f;
    char *line = NULL, *obuf, *ibuf;
    size_t cap = 0;
    ssize_t len;
    int k, selected = 0, err = 0;

    if ((o = fopen(out, append ? "a" : "w")) == NULL)
        return 0;
    obuf = ckalloc(1 << 20);
    ibuf = ckalloc(1 << 20);
    setvbuf(o, obuf, _IOFBF, 1 << 20);
    for (k = 0; k < nin; k++) {
        if ((f = fopen(in[k], "r")) == NULL) {
            fprintf(stderr, "grep: %s: %s\n", in[k], strerror(errno));
            err = 1;
            continue;
        }
        setvbuf(f, ibuf, _IOFBF, 1 << 20);
        while ((len = getline(&line, &cap, f)) != -1) {
            if (line[0] == '#' && strstr(line, pat) != NULL)
                continue;
            selected = 1;
            fwrite(line, 1, (size_t)len, o);
            if (line[len-1] != '\n')
                fputc('\n', o);
        }
        fclose(f);
    }
    free(line);
    if (fclose(o) != 0)
        err = 1;
    free(obuf);
    free(ibuf);
    return !err && selected;
}

static int ns_cat(char **in, int nin, const char *out) {
    FILE *o;
    int k, ok = 1;

    if ((o = fopen(out, "w")) == NULL)
        return 0;
    for (k = 0; k < nin; k++)
        if (!append_file(in[k], o)) {
            fprintf(stderr, "cat: %s: %s\n", in[k], strerror(errno));
            ok = 0;
        }
    if (fclose(o) != 0)
        ok = 0;
    return ok;
}

static int ns_touch(char **f, int n) {
    int k, fd, ok = 1;

    for (k = 0; k < n; k++) {
        if ((fd = open(f[k], O_WRONLY | O_CREAT, 0666)) < 0 || futimens(fd, NULL) != 0) {
            fprintf(stderr, "touch: cannot touch '%s': %s\n", f[k], strerror(errno));
            ok = 0;
        }
        if (fd >= 0)
            close(fd);
    }
    return ok;
}

static int no_options(char **w, int from, int to) {
    for (; from < to; from++)
        if (w[from][0] == '-')
            return 0;
    return 1;
}

/* echo "TEXT" > FILE or >> FILE, when TEXT holds nothing the shell or dash's echo
*  would interpret ($, `, \, "): writes TEXT and a newline, as the shell would. */
static int ns_echo(const char *cmd) {
    const char *q, *e, *t;
    FILE *o;
    int append, ok;

    if (strncmp(cmd, "echo \"", 6) != 0)
        return -1;
    q = cmd + 6;
    if ((e = strchr(q, '"')) == NULL)
        return -1;
    for (t = q; t < e; t++)
        if (*t == '$' || *t == '`' || *t == '\\' || *t == '\n')
            return -1;
    t = e + 1;
    if (strncmp(t, " >> ", 4) == 0) {
        append = 1;
        t += 4;
    } else if (strncmp(t, " > ", 3) == 0) {
        append = 0;
        t += 3;
    } else
        return -1;
    if (*t == '\0' || strpbrk(t, " \t\"'\\$`*?[]|;&<>(){}~\n") != NULL)
        return -1;
    if ((o = fopen(t, append ? "a" : "w")) == NULL)
        return 0;
    ok = fwrite(q, 1, (size_t)(e - q), o) == (size_t)(e - q) && fputc('\n', o) != EOF;
    if (fclose(o) != 0)
        ok = 0;
    return ok;
}

/* rm DIR/PREFIX*: removes every entry of DIR whose name starts with PREFIX (and,
*  like the shell's glob, not a hidden one unless PREFIX starts with '.'). With no
*  match the shell would pass the word through and rm would fail; so does this. */
static int ns_rm_glob(const char *cmd) {
    char *path, *slash, *dir, *pre;
    DIR *d;
    struct dirent *ent;
    size_t n;
    int matched = 0, ok = 1;

    if (strncmp(cmd, "rm ", 3) != 0)
        return -1;
    path = copy_string(cmd + 3);
    n = strlen(path);
    if (n < 2 || path[n-1] != '*' || strpbrk(path, " \t\"'\\$`?[]|;&<>(){}~\n") != NULL
            || strchr(path, '*') != path + n - 1 || path[0] == '-') {
        free(path);
        return -1;
    }
    path[n-1] = '\0';
    if ((slash = strrchr(path, '/')) != NULL) {
        *slash = '\0';
        dir = (slash == path) ? "/" : path;
        pre = slash + 1;
    } else {
        dir = ".";
        pre = path;
    }
    if ((d = opendir(dir)) != NULL) {
        while ((ent = readdir(d)) != NULL) {
            if (strncmp(ent->d_name, pre, strlen(pre)) != 0)
                continue;
            if (ent->d_name[0] == '.' && pre[0] != '.')
                continue;
            char *f = ckalloc(strlen(dir) + strlen(ent->d_name) + 2);
            sprintf(f, "%s/%s", dir, ent->d_name);
            matched = 1;
            if (unlink(f) != 0) {
                fprintf(stderr, "rm: cannot remove '%s': %s\n", f, strerror(errno));
                ok = 0;
            }
            free(f);
        }
        closedir(d);
    }
    if (!matched) {
        fprintf(stderr, "rm: cannot remove '%s*': No such file or directory\n", cmd[3] ? path : "");
        ok = 0;
    }
    free(path);
    return ok;
}

static int native_shell(const char *cmd) {
    char *copy, *w[NS_MAXW + 1], *out = NULL, *tok, *save = NULL;
    int nw = 0, append = 0, r = -1, k;

    if ((r = ns_echo(cmd)) >= 0 || (r = ns_rm_glob(cmd)) >= 0)
        return r;
    if (strpbrk(cmd, "\"'\\$`*?[]|;&<(){}~\n") != NULL)
        return -1;
    copy = copy_string(cmd);
    for (tok = strtok_r(copy, " \t", &save); tok != NULL; tok = strtok_r(NULL, " \t", &save)) {
        if (nw == NS_MAXW || tok[0] == '#') {
            free(copy);
            return -1;
        }
        w[nw++] = tok;
    }
    for (k = 0; k < nw; k++)               // at most one redirection: the last two words
        if (strchr(w[k], '>') != NULL) {
            if (k != nw - 2 || (strcmp(w[k], ">") != 0 && strcmp(w[k], ">>") != 0)) {
                free(copy);
                return -1;
            }
            append = (w[k][1] == '>');
            out = w[k+1];
            nw = k;
            break;
        }
    if (nw == 0) {
        free(copy);
        return -1;
    }
    w[nw] = NULL;
    if (strcmp(w[0], "mv") == 0) {
        if (nw == 3 && out == NULL && no_options(w, 1, nw)) {
            r = (rename(w[1], w[2]) == 0);
            if (!r && errno == EXDEV)
                r = copy_file(w[1], w[2]) && unlink(w[1]) == 0;
            if (!r)
                fprintf(stderr, "mv: cannot move '%s' to '%s': %s\n", w[1], w[2], strerror(errno));
        }
    } else if (strcmp(w[0], "rm") == 0) {
        if (nw >= 3 && out == NULL && strcmp(w[1], "-f") == 0 && no_options(w, 2, nw)) {
            for (k = 2; k < nw; k++)
                (void)unlink(w[k]);   // rm -f: a missing file is not an error
            r = 1;
        }
    } else if (strcmp(w[0], "touch") == 0) {
        if (nw >= 2 && out == NULL && no_options(w, 1, nw))
            r = ns_touch(w + 1, nw - 1);
    } else if (strcmp(w[0], "cat") == 0) {
        if (nw >= 2 && out != NULL && !append && no_options(w, 1, nw))
            r = ns_cat(w + 1, nw - 1, out);
    } else if (strcmp(w[0], "grep") == 0) {
        int i = 1, h = 0;
        if (i < nw && strcmp(w[i], "-v") == 0) {
            i++;
            if (i < nw && strcmp(w[i], "-h") == 0) {
                h = 1;
                i++;
            }
            if (out != NULL && i + 1 < nw && (strcmp(w[i], "eof") == 0 || strcmp(w[i], "maf") == 0)
                    && (h || nw - i - 1 == 1) && no_options(w, i + 1, nw))
                r = ns_grep_v(w[i], w + i + 1, nw - i - 1, out, append);
        }
    } else
        r = ns_spawn(w, out, append);
    free(copy);
    return r;
}
