#!/bin/bash
PD=$1
SAMPLE=$2
REF=$3
WINDOWSIZE=$4

if [ "$PD" == "poly" ]; then
   POLYDIV="-p"
elif [ "$PD" == "div" ]; then
   POLYDIV="-d"
else
   echo "Invalid site type ${PD}"
   exit 1
fi

listPolyDivSites ${POLYDIV} -n ${REF} ${SAMPLE} | nonOverlappingWindows -n -w ${WINDOWSIZE}
