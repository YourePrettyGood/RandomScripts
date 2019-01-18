#!/usr/bin/env perl

use warnings;
use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

######################################################################
# fakeHaplotype.pl                                                   #
# Usage:                                                             #
#  fakeHaplotype.pl [-i input FASTA] [-o output FASTA] [-s PRNG seed]#
#                                                                    #
# Arguments:                                                         #
#  -i,--input_FASTA         FASTA of the degenerated data to split   #
#                           (optional, defaults to STDIN)            #
#  -o,--output_FASTA        Output FASTA of the random haplotype     #
#                           (optional, defaults to STDOUT)           #
#  -s,--prng_seed           Seed for the PRNG used for splitting     #
#                           (optional, defaults to 42)               #
#                                                                    #
#  --alt_haplotype,-f       Output the other haplotype of the pair   #
#                           (flag, make sure to match the PRNG seed) #
# Description:                                                       #
#  fakeHaplotype.pl produces a fake haplotype from a degenerated     #
#  FASTA based on randomly splitting alleles at each heterozygous    #
#  site.                                                             #
######################################################################

my $SCRIPTNAME = "fakeHaplotype.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

fakeHaplotype.pl - Extract a fake haplotype from a degenerated FASTA

=head1 SYNOPSIS

fakeHaplotype.pl [options]

 Options:
  --help,-h,-?         Display this help documentation
  --input_FASTA,-i     Input FASTA of degenerate sequences to split
                       (optional, default: STDIN)
  --output_FASTA,-o    Output FASTA of one randomly-split haplotype
                       (optional, default: STDOUT)
  --prng_seed,-s       Seed to use for PRNG for splitting het sites
                       (optional, default: 42)
  --alt_haplotype,-a   Output the other haplotype of the pair
                       (flag, make sure to match the PRNG seed)
  --version,-v         Output version string

=head1 DESCRIPTION

fakeHaplotype.pl produces a fake haplotype from a degenerated
FASTA based on randomly splitting alleles at each heterozygous
site.

=cut

sub splitHaplotype($$) {
   my $line = shift @_;
   my $strand = shift @_;
   
   #Build up an array of hashes for splitting the heterozygous sites:
   my @degenerate_bases = ();
   push @degenerate_bases, {'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T', 'R'=>'A', 'Y'=>'C', 'S'=>'C', 'W'=>'A', 'K'=>'G', 'M'=>'A', 'B'=>'B', 'D'=>'D', 'H'=>'H', 'V'=>'V', 'N'=>'N', '-'=>'-', '.'=>'.'};
   push @degenerate_bases, {'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T', 'R'=>'G', 'Y'=>'T', 'S'=>'G', 'W'=>'T', 'K'=>'T', 'M'=>'C', 'B'=>'B', 'D'=>'D', 'H'=>'H', 'V'=>'V', 'N'=>'N', '-'=>'-', '.'=>'.'};
   
   #Iterate over the bases in the sequence, and split randomly at heterozygous
   # sites (do not change homozygous sites or ambiguous sites):
   my @bases = split //, $line;
   my $seq_length = scalar(@bases);
   for (my $i = 0; $i < $seq_length; $i++) {
      if ($bases[$i] !~ /[ACGTBDHVN.-]/) {
         if (int(rand(2)) == $strand) { #Choose the allele randomly
            $bases[$i] = $degenerate_bases[0]{$bases[$i]};
         } else { #Choose the other allele
            $bases[$i] = $degenerate_bases[1]{$bases[$i]};
         }
      }
   }
   return join("", @bases);
}

#Parse options:
my $help = 0;
my $man = 0;
my $in_path = '';
my $out_path = '';
my $prng_seed = 42;
my $which_haplotype = 0;
my $dispversion = 0;
GetOptions('input_FASTA|i=s' => \$in_path, 'output_FASTA|o=s' => \$out_path, 'prng_seed|s=i' => \$prng_seed, 'alt_haplotype|a' => \$which_haplotype, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

srand($prng_seed);

#Open input and output files:
my ($in, $out);
if ($in_path =~ /\.gz$/) {
   $in = new IO::Uncompress::Gunzip $in_path or die "Failed to open Gzipped input FASTA file $in_path due to error $GunzipError\n";
} elsif ($in_path eq '') {
   open $in, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input FASTA due to error $!\n";
} else {
   open $in, "<", $in_path or die "Failed to open input FASTA file due to error $!\n";
}
if ($out_path =~ /\.gz$/) {
   $out = new IO::Compress::Gzip $out_path or die "Failed to create Gzipped output FASTA file $out_path due to error $GzipError\n";
} elsif ($out_path eq '') {
   open $out, ">&", \*STDOUT or die "Failed to duplicate STDOUT file handle for output FASTA due to error $!\n";
} else {
   open $out, ">", $out_path or die "Failed to open output FASTA file due to error $!\n";
}

#Iterate over the scaffolds:
#Assumes FASTA is not line-wrapped
my $scaffold_name = ""; #Keep track of the scaffold name
while (my $line = <$in>) {
   if ($line =~ />/) {
      print $out $line;
      next;
   } else {
      chomp $line;
      #Output the random haplotype:
      my $haplotype = splitHaplotype($line, $which_haplotype);
      print $out $haplotype, "\n";
   }
}

#Close input and output files:
close $in;
close $out;

exit 0;
