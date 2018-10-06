#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to extract sites IDed by a BED file from a FASTA file
# and simply concatenate all sites under a single header.
#This is a pre-processing script for calling heterozygosity or 
# alt-homozygosity for specific sites.
#Original usage was for calling heterozygosity at 4-fold degenerate sites,
# other non-4-fold CDS sites, short intronic sites, long intronic sites,
# and intergenic sites (and perhaps UTRs).

my $SCRIPTNAME = "extractSites.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

extractSites.pl - Extract sites as a single record from a FASTA based on a BED

=head1 SYNOPSIS

extractSites.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --input_genome,-i     Path to input genome FASTA file (default: STDIN)
  --bed_path,-b         Path to BED file indicating sites to be extracted
                        from the FASTA
  --version,-v          Output version string

=head1 DESCRIPTION

This script extracts the bases at positions given by the input BED file
from the input FASTA file, and concatenates them into a single FASTA
record output to STDOUT.

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $bed_path = "";
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'bed_file|b=s' => \$bed_path, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the genome FASTA file, or set it up to be read from STDIN:
if ($genome_path ne "STDIN") {
   unless(open(GENOME, "<", $genome_path)) {
      print STDERR "Error opening genome FASTA file.\n";
      exit 2;
   }
} else {
   open(GENOME, "<&", "STDIN"); #Duplicate the file handle for STDIN to GENOME so we can seamlessly handle piping
}

#Open the BED file:
unless(open(BED, "<", $bed_path)) {
   print STDERR "Error opening BED file.\n";
   exit 3;
}

my %sites_per_scaffold = ();
while (my $line = <BED>) {
   chomp $line;
   my ($scaffold, $start_0_based, $end_1_based) = split /\t/, $line, 3;
   push @{$sites_per_scaffold{$scaffold}}, "${start_0_based}:".($end_1_based-$start_0_based);
}

close(BED);

#Now we can iterate through the genome FASTA, and concatenate sites on each
# scaffold.
my $scaffold_name = "";
my $scaffold_sequence = "";
my $sites_sequence = "";
while (my $line = <GENOME>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the sites from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "" and exists($sites_per_scaffold{$scaffold_name})) {
         #Concatenate the sites:
         for my $range (@{$sites_per_scaffold{$scaffold_name}}) {
            my ($range_start, $range_length) = split /:/, $range, 2;
            $sites_sequence .= substr $scaffold_sequence, $range_start, $range_length;
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
close(GENOME);

#Now make sure we account for the last scaffold:
if (exists($sites_per_scaffold{$scaffold_name})) {
   #Concatenate the sites:
   for my $range (@{$sites_per_scaffold{$scaffold_name}}) {
      my ($range_start, $range_length) = split /:/, $range, 2;
      $sites_sequence .= substr $scaffold_sequence, $range_start, $range_length;
   }
}

#Now output the FASTA record of sites:
print STDOUT ">Selected_sites\n", $sites_sequence, "\n";

exit 0;
