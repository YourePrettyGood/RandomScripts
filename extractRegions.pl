#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.1 (2019/04/12) Bit of extra logging                #
################################################################

#First pass script to extract regions IDed by a BED file from a FASTA file
# as individual records with headers based on the BED line.

my $SCRIPTNAME = "extractRegions.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

extractRegions.pl - Extract regions from a FASTA based on a BED

=head1 SYNOPSIS

extractRegions.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --input_genome,-i     Path to input genome FASTA file (default: STDIN)
  --bed_path,-b         Path to BED file indicating regions to be extracted
                        from the FASTA
  --prefix,-p           Prefix for FASTA headers (off by default)
  --version,-v          Output version string

=head1 DESCRIPTION

This script extracts the regions at positions given by the input BED file
from the input FASTA file, and outputs each region as its own record in a
FASTA format output to STDOUT.

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $bed_path = "";
my $header_prefix = "";
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'bed_file|b=s' => \$bed_path, 'prefix|p=s' => \$header_prefix, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;
$header_prefix .= "_" unless $header_prefix eq "";

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the genome FASTA file, or set it up to be read from STDIN:
my $genome_fh;
if ($genome_path ne "STDIN") {
   unless(open($genome_fh, "<", $genome_path)) {
      print STDERR "Error opening genome FASTA file ${genome_path}.\n";
      exit 1;
   }
} else {
   open($genome_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $genome_fh so we can seamlessly handle piping
}

#Open the BED file:
my $bed_fh;
unless(open($bed_fh, "<", $bed_path)) {
   print STDERR "Error opening BED file ${bed_path}.\n";
   exit 2;
}

my %regions_per_scaffold = ();
while (my $line = <$bed_fh>) {
   chomp $line;
   my ($scaffold, $start_0_based, $end_1_based) = split /\t/, $line, 3;
   push @{$regions_per_scaffold{$scaffold}}, "${start_0_based}:".($end_1_based-$start_0_based);
}

close($bed_fh);

#Now we can iterate through the genome FASTA, and output regions as they arise
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <$genome_fh>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the sites from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "" and exists($regions_per_scaffold{$scaffold_name})) {
         #Output the regions on this scaffold:
         for my $range (@{$regions_per_scaffold{$scaffold_name}}) {
            my ($range_start, $range_length) = split /:/, $range, 2;
            my $region_sequence = substr $scaffold_sequence, $range_start, $range_length;
            my $region_name = "";
            $region_name = ">${header_prefix}";
            $region_name .= "${scaffold_name}:" . ($range_start+1) . ".." . ($range_start+$range_length);
            print STDOUT $region_name, "\n", $region_sequence, "\n";
         }
      }
      my $scaffold_name_line = substr $line, 1; #Get rid of the prefixed ">"
      #Only preseve the first word of the scaffold header:
      my @scaffold_words = split /\s+/, $scaffold_name_line;
      $scaffold_name = $scaffold_words[0];
      $scaffold_sequence = ""; #Clear out the old sequence
   } else { #Sequence line
      $scaffold_sequence .= $line;
   }
}
close($genome_fh);

#Now make sure we account for the last scaffold:
if (exists($regions_per_scaffold{$scaffold_name})) {
   #Concatenate the sites:
   for my $range (@{$regions_per_scaffold{$scaffold_name}}) {
      my ($range_start, $range_length) = split /:/, $range, 2;
      my $region_sequence = substr $scaffold_sequence, $range_start, $range_length;
      my $region_name = "";
      $region_name = ">${header_prefix}";
      $region_name .= "${scaffold_name}:" . ($range_start+1) . ".." . ($range_start+$range_length);
      print STDOUT $region_name, "\n", $region_sequence, "\n";
   }
}

exit 0;
