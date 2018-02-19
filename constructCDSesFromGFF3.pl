#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to construct a FASTA of CDSes from a GFF3 and
# a genome FASTA.  This uses principles from fixExons.pl for
# constructing the full CDS from an Exon Range String, and
# computes the Exon Range Strings from the GFF3.

#Revised 2018/01/15 to deal with GFF3 containing more than just
# CDS records, although still assuming local sortedness
# (i.e. sorted order within a gene)

=pod

=head1 NAME

constructCDSesFromGFF3.pl - Construct a FASTA of CDSes based on a GFF3 and genome FASTA

=head1 SYNOPSIS

constructCDSesFromGFF3.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_genome,-i      Path to input genome FASTA file (default: STDIN)
  --gff3_file,-g         Path to genome annotation GFF3 file
  --output_ers,-e        Output the "Exon Range String" for the CDS?
                         0 or 1, default is 0

=head1 DESCRIPTION

This script constructs CDSes based on an input genome FASTA file and a
corresponding genome annotation GFF3 file.  The output is in a FASTA-like
format (no wrapping).

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $gff3_path = "";
my $output_exon_range_strings = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'gff3_file|g=s' => \$gff3_path, 'output_ers|e' => \$output_exon_range_strings, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

#Open the genome FASTA file, or set it up to be read from STDIN:
if ($genome_path ne "STDIN") {
   unless(open(GENOME, "<", $genome_path)) {
      print STDERR "Error opening genome FASTA file.\n";
      exit 2;
   }
} else {
   open(GENOME, "<&", "STDIN"); #Duplicate the file handle for STDIN to GENOME so we can seamlessly handle piping
}

#Open the annotation GFF3 file, which should solely contain sorted CDS records:
unless(open(GFF, "<", $gff3_path)) {
   print STDERR "Error opening GFF3 file.\n";
   exit 3;
}

#For now, we won't validate the input files, but here are our assumptions:
#FASTA file scaffold headers match up with the names in the GFF3 file
#GFF3 contains only CDS features, which are version-sorted by scaffold,
# then numerically sorted by start position
#e.g. awk '{if ($3 == "CDS") {print;};}' [original GFF3] | sort -k1,1V -k4,4n
# > [GFF3 used as input here]
#Revised to work with GFF3 containing more than just CDS records, though
# still requiring them to be locally sorted (i.e. within gene)

#Iterate over the GFF features and build the Exon Range Strings
# for each CDS:
my %CDS_coordinates = ();
my %CDSes_per_scaffold = ();
while (my $line = <GFF>) {
   chomp $line;
   next if $line =~ /^#/; #Skip comment lines
   my ($scaffold, $set, $type, $start, $end, $score, $strand, $frame, $tag_string) = split /\t/, $line, 9;
   next unless $type eq "CDS";
   my @transcript_names = ();
   if ($tag_string =~ /Parent=(.+)(?:;|$)/i) {
      @transcript_names = split /,/, $1;
   } else {
      print STDERR "Regex to find transcript name failed for tag string: ", $tag_string, "\n";
      next;
   }
   
   #Iterate over each transcript this CDS belongs to:
   for my $transcript_name (@transcript_names) {
      #exons in the Exon Range String are separated by colons
      #However, if we're at the first exon, be sure to add the transcript
      # to the list of CDSes on that scaffold
      unless (exists($CDS_coordinates{$transcript_name})) {
         push @{$CDSes_per_scaffold{$scaffold}}, $transcript_name;
      }

      #If the strand of the transcript (and thus the CDSes) is -, we denote
      # this in the exon range string by reversing the coordinates
      # (i.e. highest number first)
      if ($strand eq "+") {
         $CDS_coordinates{$transcript_name} .= ":" if exists($CDS_coordinates{$transcript_name});
         $CDS_coordinates{$transcript_name} .= "${start}-${end}";
      } else { #Need to invert the order of the exons if on - strand:
         if (exists($CDS_coordinates{$transcript_name})) {
            $CDS_coordinates{$transcript_name} = "${end}-${start}:" . $CDS_coordinates{$transcript_name};
         } else {
            $CDS_coordinates{$transcript_name} = "${end}-${start}";
         }
      }
   }
}

close(GFF);

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/ACGTNacgtn/TGCANtgcan/; #Complement
   return $reverse_sequence;
}

#As extra output, we can output a CDS range string, which translates
# between scaffold coordinate space and CDS coordinate space.
#That is, in combination with the exon range string, this provides
# a mapping between CDS position and scaffold position.
sub computeCDSRangeString($) {
   my $exon_range_string = shift @_;
   my @exons = split /:/, $exon_range_string;

   #Define the ranges in the CDS range string by recursion:
   # start_{i} = end_{i-1}+1
   # end_{i} = end_{i-1}+length(range_{i})
   # end_{0} = 0
   # length(range_{i}) can be determined from the exon range string
   #We do need to account for - strand exons in calculating length(range_{i})
   my @CDS_ranges = ();
   my $previous_range_end = 0;
   for my $range (@exons) {
      my ($exon_start_i, $exon_end_i) = split /-/, $range, 2;
      my $range_i_length;
      $range_i_length = $exon_end_i - $exon_start_i + 1 if $exon_end_i >= $exon_start_i;
      $range_i_length = $exon_start_i - $exon_end_i + 1 if $exon_start_i > $exon_end_i;
      my $CDS_start_i = $previous_range_end + 1;
      my $CDS_end_i = $previous_range_end + $range_i_length;
      push @CDS_ranges, join("-", $CDS_start_i, $CDS_end_i);
      $previous_range_end = $CDS_end_i;
   }
   return join(":", @CDS_ranges);
}

#Now we can iterate through the genome FASTA, and build the CDSes on each
# scaffold.
#Thus the output CDSes are in genome FASTA order at the scaffold level,
# but GFF order within the scaffolds.
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <GENOME>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the CDSes from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "" and exists($CDSes_per_scaffold{$scaffold_name})) {
         #Construct the CDS based on the exon range string:
         for my $CDS (@{$CDSes_per_scaffold{$scaffold_name}}) {
            my $exon_range_string = $CDS_coordinates{$CDS};
            my $constructed_CDS = "";
            my @exons = split /:/, $exon_range_string;
            for my $range (@exons) {
               my $exon;
               my ($left, $right) = split /-/, $range;
               if ($left <= $right) { #Assumes that a single base exon would never want to be reverse complemented, if it ever occurred
               #It seems a pretty safe assumption that there are no 1 bp exons.
                  my $exon_length = $right - $left + 1;
                  $exon = substr $scaffold_sequence, $left-1, $exon_length;
               } else {
                  my $exon_length = $left - $right + 1;
                  $exon = revcomp(substr($scaffold_sequence, $right-1, $exon_length));
               }
               $constructed_CDS .= $exon;
            }
            #Output the FASTA record of the CDS to STDOUT:
            print STDOUT ">", $CDS, "\n", $constructed_CDS, "\n";
            print STDOUT $scaffold_name, "=", $exon_range_string, "\n" if $output_exon_range_strings;
            print STDOUT $scaffold_name, "=", computeCDSRangeString($exon_range_string), "\n" if $output_exon_range_strings;
         }
      }
      #Only use the part of the name before the first space:
      my @scaffold_words = split /\s+/, $line;
      $scaffold_name = substr $scaffold_words[0], 1; #Get rid of the prefixed ">"
      $scaffold_sequence = ""; #Clear out the old sequence
   } else { #Sequence line
      $scaffold_sequence .= $line;
   }
}
close(GENOME);
#Now make sure we account for the last scaffold:
if (exists($CDSes_per_scaffold{$scaffold_name})) {
   for my $CDS (@{$CDSes_per_scaffold{$scaffold_name}}) {
      my $exon_range_string = $CDS_coordinates{$CDS};
      my $constructed_CDS = "";
      my @exons = split /:/, $exon_range_string;
      for my $range (@exons) {
         my $exon;
         my ($left, $right) = split /-/, $range;
         if ($left <= $right) { #Assumes that a single base exon would never want to be reverse complemented, if it ever occurred
         #It seems a pretty safe assumption that there are no 1 bp exons.
            my $exon_length = $right - $left + 1;
            $exon = substr $scaffold_sequence, $left-1, $exon_length;
         } else {
            my $exon_length = $left - $right + 1;
            $exon = revcomp(substr($scaffold_sequence, $right-1, $exon_length));
         }
         $constructed_CDS .= $exon;
      }
      #Output the FASTA record of the CDS to STDOUT:
      print STDOUT ">", $CDS, "\n", $constructed_CDS, "\n";
      print STDOUT $scaffold_name, "=", $exon_range_string, "\n" if $output_exon_range_strings;
      print STDOUT $scaffold_name, "=", computeCDSRangeString($exon_range_string), "\n" if $output_exon_range_strings;
   }
}
