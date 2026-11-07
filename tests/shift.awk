# awk -v off=N -f shift.awk in.maf > out.maf
# Moves every MAF row N bases along its source: start += N, srcSize += 2N.
# The sequence then sits N bases from either end, so forward and reverse
# strand coordinates both move by N. Drops '#' lines and collapses
# whitespace, so off=0 only normalises.
/^#/ { next }
/^s / { printf "s %s %.0f %s %s %.0f %s\n", $2, $3 + off, $4, $5, $6 + 2 * off, $7; next }
{ $1 = $1; print }
