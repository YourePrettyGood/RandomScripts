#!/bin/awk -f
/^>/{
   skip=0;
   split(substr($0, 2), idparts, "_");
#Slightly complicated regex, but basically just covers CY*, NY*, and R9b
#Might be slightly overzealous for generality, but for our data it's fine
   if (idparts[1] == "Dyak" && idparts[2] ~ /^[RCN][Y9][0-9b]+[A-C]?[0-9]?/) {
      skip=1;
   } else {
      print;
   };
}
!/^>/{
   if (skip == 0) {
      print;
   };
}
