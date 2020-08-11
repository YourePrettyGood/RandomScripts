#!/bin/bash

READFILE=$1
SORTOPTS=$2
if [[ $1 =~ .gz$ ]]; then
   READCOMMAND="gzip -dc"
#   gzip -dc ${READFILE} | awk 'NR%4==2' | sort ${SORTOPTS} | uniq -c | sort -k1,1nr
elif [[ $1 =~ .bz2$ ]]; then
   READCOMMAND="bzip2 -dc"
else
   READCOMMAND="cat"
#   awk 'NR%4==2' ${READFILE} | sort ${SORTOPTS} | uniq -c | sort -k1,1nr
fi

${READCOMMAND} ${READFILE} | awk 'NR%4==2' | sort ${SORTOPTS} | uniq -c | sort -k1,1nr
