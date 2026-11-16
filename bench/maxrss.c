/* maxrss: run a command and write the peak resident set size of the largest process in its tree (kB).
 *   maxrss OUTFILE CMD [ARGS...]
 * Uses getrusage(RUSAGE_CHILDREN) after the command exits: ru_maxrss is the largest RSS among all waited-for
 * descendants (roast's multiz children included, since roast waits for them). This avoids `time -v` inside
 * biocontainers, whose BusyBox build reports about 4x too high. Exits with the command's status. */
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: maxrss OUTFILE CMD [ARGS...]\n"); return 2; }
    pid_t p = fork();
    if (p < 0) { perror("fork"); return 2; }
    if (p == 0) { execvp(argv[2], argv + 2); perror("execvp"); _exit(127); }
    int st;
    if (waitpid(p, &st, 0) < 0) { perror("waitpid"); return 2; }
    struct rusage ru;
    getrusage(RUSAGE_CHILDREN, &ru);
    FILE *f = fopen(argv[1], "w");
    if (f) { fprintf(f, "%ld\n", ru.ru_maxrss); fclose(f); }
    return WIFEXITED(st) ? WEXITSTATUS(st) : 128 + WTERMSIG(st);
}
