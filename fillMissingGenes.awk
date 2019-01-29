#!/bin/awk -f
BEGIN{
#Make sure we parse and output tab-separated:
   FS="\t";
   OFS="\t";
}
{
#For every mRNA line in the GFF3, we check to see if the parent exists
#If the parent does not exist, we create it and output it:
   if ($3 == "mRNA") {
#Find the Parent tag, and extract the gene ID:
      split($9, tagarr, ";");
      for (tag in tagarr) {
         if (tagarr[tag] ~ /^Parent=/) {
            split(tagarr[tag], parentarr, "=");
#If the parent's never been seen before, output it:
            if (genearr[parentarr[2]] == 0) {
               print $1, $2, "gene", $4, $5, $6, $7, $8, "ID="parentarr[2];
               genearr[parentarr[2]]++;
            };
         };
#Also keep track of the mRNA records to make sure we fill them in if missing:
         if (tagarr[tag] ~ /^ID=/) {
            split(tagarr[tag], idarr, "=");
            mrnaarr[idarr[2]]++;
         };
      };
   } else {
      split($9, tagarr, ";");
      for (tag in tagarr) {
#Gene records won't have Parent tags, but exon, CDS, start_codon, and stop_codon
# records will, so if we haven't encountered the corresponding Parent for
# an exon, CDS, etc., we should create whichever parents are missing:
         if (tagarr[tag] ~ /^Parent=/) {
            split(tagarr[tag], parentarr, "=");
            split(parentarr[2], idarr, ".");
#In the case that the entire gene record is missing, we create both the gene
# and mRNA records:
            if (genearr[idarr[1]] == 0) {
               print $1, $2, "gene", $4, $5, $6, $7, $8, "ID="idarr[1];
               print $1, $2, "mRNA", $4, $5, $6, $7, $8, "ID="parentarr[2]";Parent="idarr[1];
               genearr[idarr[1]]++;
               mrnaarr[parentarr[2]]++;
            };
#In the case that the mRNA record is missing and we already created the gene
# record, we only fill in the mRNA record:
            if (mrnaarr[parentarr[2]] == 0) {
               print $1, $2, "mRNA", $4, $5, $6, $7, $8, "ID="parentarr[2]";Parent="idarr[1];
               mrnaarr[parentarr[2]]++;
            };
         };
      };
   };
#Make sure to feed-through the original line:
   print;
}
