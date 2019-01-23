#!/bin/awk -f
/^>/{
#Could modify this so that we append something to the header to
# indicate that it has been trimmed:
   print;
}
!/^>/{
#This line does quite a lot at once:
#It trims length % 3 bases off the right end of the sequence by
# taking the leftmost prefix of length 3*floor(length/3)
# and it converts this to all uppercase:
   trimmed=toupper(substr($0, 1, 3*int(length($0)/3)));
#Now to get rid of terminal stop codons, we use a regex based on
# the universal genetic code (NCBI transl_table=1), using a
# zero-width assertion for the end of the line to ensure it
# is a terminal stop codon that is trimmed, and using sub()
# to ensure it's only trimmed once:
   sub(/T(AA|AG|GA)$/, "", trimmed);
   print trimmed;
}
