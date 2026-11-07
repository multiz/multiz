#ifndef SPECIESTREE_H
#define SPECIESTREE_H
#include "align_util.h"

void do_cmd(const char *fmt, ...);

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

/* parse treeStr, tree is allocated, assuming enough size */
int  parseSpeciesTree(char* treeStr, TreeNode* tree, int nbz, char** bz_file, int* id, char* buf, int (*operation)(TreeNodePtr, TreeNodePtr, int, int, char**));

//NameListPtr formNameList(char* str);

#endif
