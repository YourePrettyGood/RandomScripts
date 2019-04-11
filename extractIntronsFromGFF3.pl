#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.0 (2019/04/09) Initial version                     #
################################################################

#First pass script to construct a FASTA of introns from a GFF3 and
# a genome FASTA.  This uses structures from constructCDSesFromGFF3.pl
# and adds in flanking exon lengths, IDs, and intron length to
# facilitate downstream filtration for homology.

#Assumes GFF3 is locally sorted (i.e. sorted order within a gene)
#Output order is guaranteed to match two conditions:
# 1) Sortedness within gene (increasing coordinate order)
# 2) Scaffold order is as found in the genome FASTA

my $SCRIPTNAME = "extractIntronsFromGFF3.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

extractIntronsFromGFF3.pl - Construct a FASTA of introns based on a GFF3 and genome FASTA

=head1 SYNOPSIS

extractIntronsFromGFF3.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_genome,-i      Path to input genome FASTA file (default: STDIN)
  --gff3_file,-g         Path to genome annotation GFF3 file
  --version,-v           Output version string

=head1 DESCRIPTION

This script constructs a FASTA of introns from a GFF3 and a genome FASTA.
This uses structures from constructCDSesFromGFF3.pl and adds a header with
intron ID (GFF3 format), intron coordinates (FlyBase-style), intron length,
flanking exon lengths, and flanking exon IDs to facilitate downstream
filtration for homology.

Assumes GFF3 is locally sorted (i.e. sorted order within a gene)
Output order is guaranteed to match two conditions:
 1) Sortedness within gene (increasing coordinate order)
 2) Scaffold order is as found in the genome FASTA

=cut

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/AaCcGgTtRrYySsWwKkMmBbDdHhVvNn/TtGgCcAaYyRrSsWwMmKkVvHhDdBbNn/; #Complement incl. IUPAC degenerate bases
   return $reverse_sequence;
}

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $gff3_path = "";
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'gff3_file|g=s' => \$gff3_path, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
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

#Open the annotation GFF3 file:
unless(open(GFF, "<", $gff3_path)) {
   print STDERR "Error opening GFF3 file.\n";
   exit 3;
}

#For now, we won't validate the input files, but here are our assumptions:
#FASTA file scaffold headers match up with the names in the GFF3 file
#GFF3 features are version-sorted by scaffold, then numerically sorted
# by start position
#e.g. sort -k1,1V -k4,4n < [original GFF3] > [GFF3 used as input here]

my $FASTA_skip = 0; #Skip all lines after the ##FASTA line if present

#Iterate over the GFF features and store exons into our data structure:
my %exon_coordinates = ();
my %exon_IDs = ();
my %tx_per_scaffold = ();
while (my $line = <GFF>) {
   chomp $line;
   $FASTA_skip = 1 if $line =~ /^##FASTA/; #Trigger FASTA skip
   next if $FASTA_skip;
   next if $line =~ /^#/; #Skip comment lines
   my ($scaffold, $set, $type, $start, $end, $score, $strand, $frame, $tag_string) = split /\t/, $line, 9;
   next unless $type eq "exon";
   my @transcript_names = ();
   if ($tag_string =~ /Parent=(.+?)(?:;|$)/i) {
      @transcript_names = split /,/, $1;
   } else {
      print STDERR "Regex to find transcript name failed for tag string: ", $tag_string, "\n";
      next;
   }
   my $exon_ID = "";
   if ($tag_string =~ /ID=(.+?)(?:;|$)/i) {
      $exon_ID = $1;
   } else {
      print STDERR "Regex to find exon ID failed for tag string: ", $tag_string, "\n";
      next;
   }
   
   #Iterate over each transcript this exon belongs to:
   for my $transcript_name (@transcript_names) {
      #exons in the Exon Range String are separated by colons
      #However, if we're at the first exon, be sure to add the transcript
      # to the list of exons on that scaffold
      unless (exists($exon_coordinates{$transcript_name})) {
         push @{$tx_per_scaffold{$scaffold}}, $transcript_name;
      }

      #If the strand of the transcript (and thus the exons) is -, we denote
      # this in the exon range string by reversing the coordinates
      # (i.e. highest number first)
      if ($strand eq "+") {
         $exon_coordinates{$transcript_name} .= ":" if exists($exon_coordinates{$transcript_name});
         $exon_coordinates{$transcript_name} .= "${start}-${end}";
         $exon_IDs{$transcript_name} .= ":" if exists($exon_IDs{$transcript_name});
         $exon_IDs{$transcript_name} .= ${exon_ID};
      } else { #Need to invert the order of the exons if on - strand:
         if (exists($exon_coordinates{$transcript_name})) {
            $exon_coordinates{$transcript_name} = "${end}-${start}:" . $exon_coordinates{$transcript_name};
            $exon_IDs{$transcript_name} = "${exon_ID}:" . $exon_IDs{$transcript_name};
         } else {
            $exon_coordinates{$transcript_name} = "${end}-${start}";
            $exon_IDs{$transcript_name} = ${exon_ID};
         }
      }
   }
}

close(GFF);

#Now we can iterate through the genome FASTA, and extract the introns
#Thus the output introns are in genome FASTA order at the scaffold level,
# but GFF order within the scaffolds.
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <GENOME>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the introns from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "" and exists($tx_per_scaffold{$scaffold_name})) {
         #Extract introns based on the exon range string:
         for my $transcript (@{$tx_per_scaffold{$scaffold_name}}) {
            my $exon_range_string = $exon_coordinates{$transcript};
            my @exons = split /:/, $exon_range_string;
            my @IDs = split /:/, $exon_IDs{$transcript};
            my $num_exons = scalar(@exons);
            #Iterate over the introns by iterating over exons, skipping the first:
            for (my $i = 1; $i < $num_exons; $i++) {
               my ($prev_left, $prev_right) = split /-/, $exons[$i-1], 2;
               my ($cur_left, $cur_right) = split /-/, $exons[$i], 2;
               my $prev_length = abs($prev_right - $prev_left) + 1;
               my $cur_length = abs($cur_right - $cur_left) + 1;
               my $prev_ID = $IDs[$i-1];
               my $cur_ID = $IDs[$i];
               my $intron_length = abs($cur_left - $prev_right) - 1; # -1 = + 1 - 2 because the coordinates are for the flanking bases of the exons, so we skip that flanking 1 bp on either side to get intron length
               my $intron_seq = "";
               my $intron_coords = "";
               my $intron_strand = "+";
               next if $intron_length < 1; #This should never occur, at the very least we expect splice junctions, so maybe this should be < 4
               if ($prev_right < $cur_left) { #Forward strand
                  #On forward strand, prev_right is equal to 0-based lower flank of intron
                  $intron_seq = substr $scaffold_sequence, $prev_right, $intron_length;
                  $intron_coords = join("..", $prev_right+1, $cur_left-1);
               } else { #Reverse strand
                  #On reverse strand, cur_left is equal to 0-based lower flank of intron
                  $intron_seq = revcomp(substr($scaffold_sequence, $cur_left, $intron_length));
                  $intron_coords = join("..", $cur_left+1, $prev_right-1);
                  $intron_strand = "-";
               }
               #Print the FASTA header (Intron ID, coords, Intron Length, Left Exon Length, Right Exon Length, Left Exon ID, Right Exon ID):
               print STDOUT ">${transcript}.intron${i} ${scaffold_name}:${intron_coords}(${intron_strand}) ${intron_length} ${prev_length} ${cur_length} ${prev_ID} ${cur_ID}\n";
               #Now print the intron sequence:
               print STDOUT "${intron_seq}\n";
            }
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
if (exists($tx_per_scaffold{$scaffold_name})) {
   #Extract introns based on the exon range string:
   for my $transcript (@{$tx_per_scaffold{$scaffold_name}}) {
      my $exon_range_string = $exon_coordinates{$transcript};
      my @exons = split /:/, $exon_range_string;
      my @IDs = split /:/, $exon_IDs{$transcript};
      my $num_exons = scalar(@exons);
      #Iterate over the introns by iterating over exons, skipping the first:
      for (my $i = 1; $i < $num_exons; $i++) {
         my ($prev_left, $prev_right) = split /-/, $exons[$i-1], 2;
         my ($cur_left, $cur_right) = split /-/, $exons[$i], 2;
         my $prev_length = abs($prev_right - $prev_left) + 1;
         my $cur_length = abs($cur_right - $cur_left) + 1;
         my $prev_ID = $IDs[$i-1];
         my $cur_ID = $IDs[$i];
         my $intron_length = abs($cur_left - $prev_right) - 1; # -1 = + 1 - 2 because the coordinates are for the flanking bases of the exons, so we skip that flanking 1 bp on either side to get intron length
         my $intron_seq = "";
         my $intron_coords = "";
         my $intron_strand = "+";
         next if $intron_length < 1; #This should never occur, at the very least we expect splice junctions, so maybe this should be < 4
         if ($prev_right < $cur_left) { #Forward strand
            $intron_seq = substr $scaffold_sequence, $prev_right+1, $intron_length;
            $intron_coords = join("..", $prev_right+1, $cur_left-1);
         } else { #Reverse strand
            $intron_seq = revcomp(substr($scaffold_sequence, $prev_right+1, $intron_length));
            $intron_coords = join("..", $cur_left+1, $prev_right-1);
            $intron_strand = "-";
         }
         #Print the FASTA header (Intron ID, coords, Intron Length, Left Exon Length, Right Exon Length, Left Exon ID, Right Exon ID):
         print STDOUT ">${transcript}.intron${i} ${scaffold_name}:${intron_coords}(${intron_strand}) ${intron_length} ${prev_length} ${cur_length} ${prev_ID} ${cur_ID}\n";
         #Now print the intron sequence:
         print STDOUT "${intron_seq}\n";
      }
   }
}

exit 0;
