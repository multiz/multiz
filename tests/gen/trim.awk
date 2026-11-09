# trim.awk: keep the first NB blocks, the first NR rows of each, and drop columns
# that are '-' in every kept row. Usage: awk -v nb=2 -v nr=40 -f trim.awk in.maf
function flush(   i, j, keep, t, out) {
    if (nrow == 0) return
    for (j = 1; j <= len; j++) {
        keep[j] = 0
        for (i = 1; i <= nrow; i++) if (substr(txt[i], j, 1) != "-") { keep[j] = 1; break }
    }
    print aline
    for (i = 1; i <= nrow; i++) {
        out = ""
        for (j = 1; j <= len; j++) if (keep[j]) out = out substr(txt[i], j, 1)
        print pre[i] out
    }
    print ""
    nrow = 0
}
/^#/ { if (blk == 0) print; next }
/^a/ { flush(); blk++; if (blk > nb) exit; aline = $0; next }
/^s/ { if (blk > nb || nrow >= nr) next
       t = $NF; p = substr($0, 1, length($0) - length(t))
       nrow++; txt[nrow] = t; pre[nrow] = p; len = length(t); next }
END { flush() }
