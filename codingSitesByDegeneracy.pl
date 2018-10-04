#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#For debugging:
#use Data::Dumper;

################################################################
#                                                              #
################################################################

#Script to identify n-fold degenerate sites (n \in {1,2,3,4}) from
# a FASTA of CDSes.
#If the FASTA was extracted from the genome+GFF3 using
# constructCDSesFromGFF3.pl, the BED intervals output by this
# script may be sorted, merged, and then fed into CDStoGenomicIntervals.pl
# in order to obtain a BED in genomic coordinates of these sites.
#This can facilitate subsetting genome-wide Dxy or Pi produced by
# calculateDxy or calculatePolymorphism (or listPolyDivSites)
# using subsetVCFstats.pl.

=pod

=head1 NAME

codingSitesByDegeneracy.pl - Identifies n-fold degenerate sites in a FASTA of CDSes

=head1 SYNOPSIS

CDStoGenomicIntervals.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_FASTA,-i       Path to input FASTA of CDSes
                         (default: STDIN)
  --degeneracy,-f        Degree of degeneracy of site to output
                         (i.e. 1, 2, 3, or 4)
                         (1-fold aka 0-fold)
  --debug,-d             Extra debugging output

=head1 DESCRIPTION

This script outputs BED intervals (unsorted, unmerged) corresponding to
sites in each CDS that have a certain level of degeneracy (specified by
the -f argument).
For instance, if you want to identify all 4-fold degenerate synonymous
sites in a FASTA of CDSes called myCDSes.fasta, you would do:
codingSitesByDegeneracy.pl -i myCDSes.fasta -f 4 | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > my4foldsites.bed

For 0-fold degenerate sites (i.e. sites where any mutational direction would
change the amino acid encoded by the codon), specify -f 0 or -f 1.

In the exceedingly rare case you are interested in 2-fold or 3-fold degenerate
sites, the script can handle it.  Just don't expect very many sites.

=cut

my $help = 0;
my $man = 0;
my $cds_path = "STDIN";
my $degeneracy = 4;
my $debug = 0;
GetOptions('input_FASTA|i=s' => \$cds_path, 'degeneracy|f=i' => \$degeneracy, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

#Catch the terminology case (0-fold should really be called 1-fold):
$degeneracy = 1 if $degeneracy == 0;
#Catch the invalid input case:
print STDERR "You input an invalid degree of degeneracy (${degeneracy}), using 4 instead\n" if $degeneracy < 0 or $degeneracy > 4;
$degeneracy = 4 if $degeneracy < 0 or $degeneracy > 4;

#Some basic constants:
my $NUM_BASES = 4;
my $CODON_LENGTH = 3;

#Determine degeneracy for every possible site in every possible codon:
my @translate = ("K", "Q", "E", "*", "T", "P", "A", "S", "R", "R", "G", "*", "I", "L", "V", "L", "N", "H", "D", "Y", "T", "P", "A", "S", "S", "R", "G", "C", "I", "L", "V", "F", "K", "Q", "E", "*", "T", "P", "A", "S", "R", "R", "G", "W", "M", "L", "V", "L", "N", "H", "D", "Y", "T", "P", "A", "S", "S", "R", "G", "C", "I", "L", "V", "F");
#We encode codons as integers from 0-63 here, in least significant bit order
# so AAA=000000=0, TTT=111111=63, AAC=000001=1, CAA=010000=16, AAG=000010=2
#In order to identify all possible alternate codons at a given position, we
# simply toggle one or both bits of the appropriate position.
my @degen_positions = ();
for (my $i = 0; $i < $NUM_BASES**$CODON_LENGTH; $i++) {
   $degen_positions[$i] = [];
   #Iterate through codon positions:
   for (my $j = 0; $j < $CODON_LENGTH; $j++) {
      #Minimum degree of degeneracy is necessarily 1:
      my $degen_count = 1;
      #Iterate through alternate codons:
      for (my $k = 1; $k < $NUM_BASES; $k++) {
         $degen_count++ if $translate[$i] eq $translate[$i ^ ($k * $NUM_BASES ** $j)];
      }
      push @{$degen_positions[$i]}, $degen_count;
   }
}

#Debug @degen_positions:
#print STDERR Dumper(\@degen_positions) if $debug;
#print STDERR scalar(@degen_positions) if $debug;
#exit 1 if $debug > 1;

#Note: This accounts for degeneracy of stop codons, though they should be
# excluded from CDS FASTAs, typically

#We need a function to convert codon character strings to our integer
# representation.
#Any codon containing non-ACGT characters is set to 64, and filtered
# out of any further processing.
sub codon_to_int($) {
   my $codon = shift @_;
   return 64 if $codon !~ /[ACGTacgt]+/;
   my %nuc_to_int = ('A' => 0, 'C' => 1, 'G' => 2, 'T' => 3);
   my @nucs = split //, $codon;
   print STDERR "Invalid codon of length ", length($codon), "\n" unless length($codon) == 3;
   return 65 unless length($codon) == 3;
   my $codon_int = 0;
   for (my $i = 2; $i >= 0; $i--) {
      $codon_int += $nuc_to_int{uc($nucs[$i])} * 4 ** $i;
   }
   return $codon_int;
}

#Open the CDS FASTA file, or set it up to be read from STDIN:
my $cdsfh;
if ($cds_path ne "STDIN") {
   unless(open($cdsfh, "<", $cds_path)) {
      print STDERR "Error opening CDS FASTA file ${cds_path}.\n";
      exit 2;
   }
} else {
   open($cdsfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $cdsfh so we can seamlessly handle piping
}

#Read the CDSes into a hash for easy access and iteration:
print STDERR "Reading CDS FASTA ${cds_path}\n";
my %CDSes = ();
my ($header, $seq) = ("", "");
while (my $line = <$cdsfh>) {
   chomp $line;
   if ($line =~ /^>/) {
      if ($header ne "") {
         $CDSes{$header} = $seq;
         $seq = "";
      }
      #We ignore anything after a whitespace character, as those headers get annoying...
      my @headerparts = split /\s+/, $line;
      $header = substr($headerparts[0], 1);
   } else {
      $seq .= $line;
   }
}
#Make sure to store the last sequence:
$CDSes{$header} = $seq;

close($cdsfh);

#Now iterate over each CDS to identify degenerate sites:
print STDERR "Identifying ${degeneracy}-fold degenerate sites in these CDSes\n";
for my $CDS (keys %CDSes) {
   my $CDS_seq = $CDSes{$CDS};
   my $CDS_length = length($CDS_seq);
   print STDERR "CDS length for ${CDS} not a multiple of 3, it is ${CDS_length}\n" unless $CDS_length % 3 == 0;
   next unless $CDS_length % 3 == 0;
   for (my $i = 0; $i < $CDS_length; $i += 3) {
      my $codon = codon_to_int(substr($CDS_seq, $i, 3));
      print STDERR "Unable to process codon in ${CDS} starting at position " . $i+1 . ", as it contains a non-ACGT base\n" if $codon < 0 or $codon > 63;
      next if $codon < 0 or $codon > 63;
      for (my $j = 0; $j < scalar(@{$degen_positions[$codon]}); $j++) {
         #Debugging:
         #print STDERR "Degeneracy ${degeneracy}\n" if $debug > 1;
         #print STDERR "Codon ${codon}\n" if $debug > 1;
         #print STDERR "Base in codon $j\n" if $debug > 1;
         #print STDERR Dumper($degen_positions[$codon]) if $debug > 1;
         #Print the BED interval for this position:
         print join("\t", $CDS, $i+$j, $i+$j+1), "\n" if $degen_positions[$codon][$j] == $degeneracy;
      }
   }
}

exit 0;
