#ifndef SPECIESTREE_H
#define SPECIESTREE_H
#include "align_util.h"

void do_cmd(const char *fmt, ...);

/* roast and tba build shell-style commands from user paths, and do_cmd hands any it
*  cannot run natively to /bin/sh unquoted. These refuse, up front, what the shell would
*  split or interpret, so a user path never reaches /bin/sh as more than one plain word.
*  check_shell_word: a path or name (no white space, quotes, $ ` \ globs, | ; & < > ( ) { } ~,
*  nor a leading -). check_shell_text: text echoed into the header (no " $ ` \ or newline). */
void check_shell_word(const char *what, const char *s);

/* make_work_dir: roast's and tba's private temporary directory, DIR/_TAG_PID_XXXXXX (mkdtemp), returned as a
*  prefix ending in '/' for the run's file names; with execute 0 (a dry run) nothing is created. Two runs that
*  share a T= (NFS, or containers with the same small pids) used to share file names, overwriting and deleting
*  each other's files. remove_work_dir removes the (by then empty) directory. */
char *make_work_dir(const char *dir, const char *tag);
void remove_work_dir(void);
void check_shell_text(const char *what, const char *s);

/* The file steps roast and tba run between merges, done natively instead of through
*  system(). Each prints the command it stands for when verbose is set, does nothing
*  when execute is 0, and on failure stops only if force is 0, exactly as do_cmd does
*  for the shell command it replaces.
*/
void native_mv(const char *from, const char *to);                  // mv from to
void native_rm2(const char *a, const char *b);                      // rm -f a b
void native_cp(const char *from, const char *to);                   // cp from to
void native_grep_v_append(const char *pat, const char *in1, const char *in2, const char *out);
    // grep -v pat in1 >> out, or grep -v -h pat in1 in2 >> out when in2 is not NULL
void native_run(const char *words, const char *const *args, const char *out, int append);
    /* words args... > out (append: >> out), started directly rather than through
    *  /bin/sh. words is a program name plus options, split on spaces as the shell
    *  would; args follow it. */

/* parse treeStr, tree is allocated, assuming enough size */
int  parseSpeciesTree(char* treeStr, TreeNode* tree, int nbz, char** bz_file, int* id, int (*operation)(TreeNodePtr, TreeNodePtr, int, int, char**));

//NameListPtr formNameList(char* str);

#endif
