#!/bin/awk -f
#We assume that the input is not line-wrapped
BEGIN{
   i = 0;
   minlen = 999999999; #You'll probably never find this in an alignment
}
/^>/{
   header[i] = $0;
}
!/^>/{
   #We store the minimum length
   minlen = length($0) < minlen ? length($0) : minlen;
   sequence[i] = toupper($0);
   i += 1;
}
END{
   for (j = 0; j < i; j += 1) {
      print header[j];
      print substr(sequence[j], 1, minlen);
   }
}
