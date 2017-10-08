#!/bin/bash

READFILE=$1
if [[ $1 =~ .gz$ ]]; then
   zcat ${READFILE} | awk 'NR%4==2' | sort | uniq -c | sort -k1,1nr
else
   awk 'NR%4==2' ${READFILE} | sort | uniq -c | sort -k1,1nr
fi
