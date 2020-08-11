#!/bin/awk -f
#Assumptions:
#Input GFF3 is *not* missing any mRNA records, and mRNA records
# have a Parent tag in column 9
#This avoids needing to infer the Parent ID
#We also don't output inferred parent gene records inline, but
# instead output them at the end, to make sure we have the
# correct interval (since the first mRNA may not be the longest)
BEGIN{
#Make sure we parse and output tab-separated:
   FS="\t";
   OFS="\t";
}
#Feed through any header lines:
/^#/{
   print;
}
!/^#/{
   split($9, tags, ";");
   for (i in tags) {
      split(tags[i], elems, "=");
      if (elems[1] == "ID") {
         id=elems[2];
      } else if (elems[1] == "Parent") {
         parent=elems[2];
      }
   }
   if ($3 == "gene") {
   #Store any gene records for update and output at end:
      if (length(genelefts[id]) == 0) {
         genelefts[id]=$1"\t"$2"\t"$3;
         genestarts[id]=$4;
         geneends[id]=$5;
         genescores[id]=$6;
         genestrands[id]=$7;
         generights[id]=$8"\t"$9;
      } else {
         print "Duplicated gene record found!" > "/dev/stderr";
         print genelefts[id], genestarts[id], geneends[id], genescores[id], genestrands[id], generights[id] > "/dev/stderr";
         print $0 > "/dev/stderr";
      }
   } else if ($3 == "mRNA") {
   #Check if mRNA record's parent exists in gene set, create if not:
      if (length(genelefts[parent]) == 0) {
         genelefts[parent]=$1"\t"$2"\tgene";
         genestarts[parent]=$4;
         geneends[parent]=$5;
         genescores[parent]=$6;
         genestrands[parent]=$7;
         generights[parent]=$8"\tID="parent;
      } else {
      #Parent gene exists, so just update starts and ends if necessary:
         if ($4 < genestarts[parent]) {
            genestarts[parent]=$4;
         }
         if ($5 > geneends[parent]) {
            geneends[parent]=$5;
         }
      }
      print $0;
   } else {
   #Feed through any records other than mRNA or gene:
      print $0;
   }
}
END{
   for (i in genelefts) {
      print genelefts[i], genestarts[i], geneends[i], genescores[i], genestrands[i], generights[i];
   }
}
#For every mRNA line in the GFF3, we check to see if the parent exists
#If the parent does not exist, we create it and output it:
#   if ($3 == "mRNA") {
#Find the Parent tag, and extract the gene ID:
#      split($9, tagarr, ";");
#      for (tag in tagarr) {
#         if (tagarr[tag] ~ /^Parent=/) {
#            split(tagarr[tag], parentarr, "=");
#If the parent's never been seen before, output it:
#            if (genearr[parentarr[2]] == 0) {
#               print $1, $2, "gene", $4, $5, $6, $7, $8, "ID="parentarr[2];
#               genearr[parentarr[2]]++;
#            };
#         };
#Also keep track of the mRNA records to make sure we fill them in if missing:
#         if (tagarr[tag] ~ /^ID=/) {
#            split(tagarr[tag], idarr, "=");
#            mrnaarr[idarr[2]]++;
#         };
#      };
#   } else {
#      split($9, tagarr, ";");
#      for (tag in tagarr) {
#Gene records won't have Parent tags, but exon, CDS, start_codon, and stop_codon
# records will, so if we haven't encountered the corresponding Parent for
# an exon, CDS, etc., we should create whichever parents are missing:
#         if (tagarr[tag] ~ /^Parent=/) {
#            split(tagarr[tag], parentarr, "=");
#            split(parentarr[2], idarr, ".");
#In the case that the entire gene record is missing, we create both the gene
# and mRNA records:
#            if (genearr[idarr[1]] == 0) {
#               print $1, $2, "gene", $4, $5, $6, $7, $8, "ID="idarr[1];
#               print $1, $2, "mRNA", $4, $5, $6, $7, $8, "ID="parentarr[2]";Parent="idarr[1];
#               genearr[idarr[1]]++;
#               mrnaarr[parentarr[2]]++;
#            };
#In the case that the mRNA record is missing and we already created the gene
# record, we only fill in the mRNA record:
#            if (mrnaarr[parentarr[2]] == 0) {
#               print $1, $2, "mRNA", $4, $5, $6, $7, $8, "ID="parentarr[2]";Parent="idarr[1];
#               mrnaarr[parentarr[2]]++;
#            };
#         };
#      };
#   };
#Make sure to feed-through the original line:
#   print;
#}
