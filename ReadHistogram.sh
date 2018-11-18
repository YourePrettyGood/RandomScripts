#!/bin/bash

READFILE=$1
SORTOPTS=$2
if [[ $1 =~ .gz$ ]]; then
   zcat ${READFILE} | awk 'NR%4==2' | sort ${SORTOPTS} | uniq -c | sort -k1,1nr
else
   awk 'NR%4==2' ${READFILE} | sort ${SORTOPTS} | uniq -c | sort -k1,1nr
fi
