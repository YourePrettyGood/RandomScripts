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

listPolyDivSites ${POLYDIV} -n ${REF} ${SAMPLE} | awk 'BEGIN{OFS="\t";}{if($4==0){stat[$1]+=$3;count[$1]+=1;}}END{for (scaf in stat){if (count[scaf] > 0){print scaf,stat[scaf]/count[scaf];}else{print scaf,"NA";}}}' | sort -k1,1V
