#!/bin/bash
PD=$1
SAMPLE=$2
REF=$3

if [ "$PD" == "poly" ]; then
   POLYDIV="-p"
elif [ "$PD" == "div" ]; then
   POLYDIV="-d"
else
   echo "Invalid site type ${PD}"
   exit 1
fi

listPolyDivSites ${POLYDIV} -n ${REF} ${SAMPLE} | awk 'BEGIN{OFS="\t";}{if($4==0){stat+=$3;count+=1;}}END{if (count > 0) {print stat/count;}else{print "NA";}}'
