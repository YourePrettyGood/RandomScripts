#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to translate BED intervals in CDS space to
# genomic space based on a GFF3
#Assumes sorted order of records within a gene in a GFF3

=pod

=head1 NAME

CDStoGenomicIntervals.pl - Convert CDS-space BED intervals to genomic-space using a GFF3

=head1 SYNOPSIS

CDStoGenomicIntervals.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_BED,-i         Path to input BED file in CDS-space
                         (default: STDIN)
  --gff3_file,-g         Path to genome annotation GFF3 file
  --debug,-d             Extra debugging output

=head1 DESCRIPTION

This script converts BED intervals in CDS space to BED intervals in
genomic space using a GFF3.

=cut

my $help = 0;
my $man = 0;
my $bed_path = "STDIN";
my $gff3_path = "";
my $debug = 0;
GetOptions('input_BED|i=s' => \$bed_path, 'gff3_file|g=s' => \$gff3_path, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

#Open the CDS-space BED file, or set it up to be read from STDIN:
my $bedfh;
if ($bed_path ne "STDIN") {
   unless(open($bedfh, "<", $bed_path)) {
      print STDERR "Error opening CDS-space BED file ${bed_path}.\n";
      exit 2;
   }
} else {
   open($bedfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $bedfh so we can seamlessly handle piping
}

#Open the annotation GFF3 file, which should solely contain sorted CDS records:
my $gff_fh;
unless(open($gff_fh, "<", $gff3_path)) {
   print STDERR "Error opening GFF3 file ${gff3_path}.\n";
   exit 3;
}

#For now, we won't validate the input files, but here are our assumptions:
#BED file scaffold IDs match up with transcript IDs in the GFF3 file
#GFF3 features are version-sorted by scaffold, then numerically sorted
# by start position
#e.g. sort -k1,1V -k4,4n < [original GFF3] > [GFF3 used as input here]

my $FASTA_skip = 0; #Skip all lines after the ##FASTA line if present

#Iterate over the GFF features and build the Exon Range Strings
# for each CDS:
print STDERR "Processing ${gff3_path} and building genome-space exon interval arrays\n";
my %ers = ();
my %transcript_scaffold_map = ();
while (my $line = <$gff_fh>) {
   chomp $line;
   $FASTA_skip = 1 if $line =~ /^##FASTA/; #Trigger FASTA skip
   next if $FASTA_skip;
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
      $transcript_scaffold_map{$transcript_name} = {"scaffold" => $scaffold, "strand" => $strand} unless exists($transcript_scaffold_map{$transcript_name});
      #Initialize the array of exons in %ers for this transcript if it's new:
      $ers{$transcript_name} = [] unless exists($ers{$transcript_name});
      #Add the exon to %ers:
      #If the strand of the transcript (and thus the CDSes) is -, we
      # reverse the coordinate order to keep alignment with CDS-space
      # intervals (i.e. highest number first)
      if ($strand eq "+") {
         my $exon = {"start" => $start, "end" => $end};
         push @{$ers{$transcript_name}}, $exon;
      } else { #Need to invert the order of the exons if on - strand:
         my $exon = {"start" => $end, "end" => $start};
         unshift @{$ers{$transcript_name}}, $exon;
      }
   }
}

close($gff_fh);

#Now create a corresponding CDS range array, which will correspond to
# CDS-space rather than genomic-space:
print STDERR "Generating CDS-space intervals for transcripts\n";
my %cds = ();
for my $transcript (keys %ers) {
   my $strand = $transcript_scaffold_map{$transcript}{'strand'};
   $cds{$transcript} = [] unless exists($cds{$transcript});
   my $previous_range_end = 0;
   for my $exon (@{$ers{$transcript}}) {
      my $range_length;
      $range_length = $exon->{'end'} - $exon->{'start'} + 1 if $strand eq "+";
      $range_length = $exon->{'start'} - $exon->{'end'} + 1 unless $strand eq "+";
      my $CDS_start = $previous_range_end + 1;
      my $CDS_end = $previous_range_end + $range_length;
      $previous_range_end = $CDS_end;
      push @{$cds{$transcript}}, {"start" => $CDS_start, "end" => $CDS_end};
   }
}

#Now read in the BED intervals and convert them appropriately:
print STDERR "Reading ${bed_path} and converting intervals to genome-space\n";
while (my $line = <$bedfh>) {
   chomp $line;
   my ($transcript, $BEDstart, $BEDend) = split /\t/, $line, 3;
   my %query = ("start" => $BEDstart+1, "end" => $BEDend);
   print STDERR "Contents of %query: ", $query{'start'}, " ", $query{'end'}, "\n" if $debug > 1;
   unless (exists($transcript_scaffold_map{$transcript}) and exists($ers{$transcript})) {
      print STDERR "Skipping transcript, unable to find a GFF3 entry for transcript ${transcript}\n";
      next;
   }
   my $scaffold = $transcript_scaffold_map{$transcript}{'scaffold'};
   my $strand = $transcript_scaffold_map{$transcript}{'strand'};
   my ($start_exon_index, $end_exon_index) = (0, 0);
   #Identify bounding exons:
   print STDERR "SEI: ", $start_exon_index, " EEI: ", $end_exon_index, " Transcript: ", $transcript, " # exons: ", scalar(@{$cds{$transcript}}), "/", scalar(@{$ers{$transcript}}), "\n" if $debug > 1;
   $start_exon_index++ until $query{'start'} >= $cds{$transcript}[$start_exon_index]{'start'} and $query{'start'} <= $cds{$transcript}[$start_exon_index]{'end'};
   $end_exon_index++ until $query{'end'} >= $cds{$transcript}[$end_exon_index]{'start'} and $query{'end'} <= $cds{$transcript}[$end_exon_index]{'end'};
   print STDERR "SEI: ", $start_exon_index, " EEI: ", $end_exon_index, " Transcript: ", $transcript, " # exons: ", scalar(@{$cds{$transcript}}), "/", scalar(@{$ers{$transcript}}), "\n" if $debug > 1;
   #Make sure to account for intervals spanning multiple exons:
   while ($start_exon_index <= $end_exon_index) {
      #Calculate offset of interval start from start of exon:
      my $start_dist = $query{'start'} - $cds{$transcript}[$start_exon_index]{'start'};
      $start_dist = -$start_dist unless $strand eq "+";
      #Calculate genome-space coordinates of interval start:
      my $genomic_start = $ers{$transcript}[$start_exon_index]{'start'} + $start_dist;
      #Move the query start so that we don't duplicate output intervals:
      $query{'start'} = $cds{$transcript}[$start_exon_index]{'end'}+1;
      #Calculate positive offset of interval end from end of exon:
      my $end_dist = $query{'end'} < $cds{$transcript}[$start_exon_index]{'end'} ? $cds{$transcript}[$start_exon_index]{'end'} - $query{'end'} : 0;
      $end_dist = -$end_dist unless $strand eq "+";
      #Calculate genome-space coordinates of interval end:
      my $genomic_end = $ers{$transcript}[$start_exon_index]{'end'} - $end_dist;
      #Now print the interval:
      print join("\t", $scaffold, $genomic_start-1, $genomic_end), "\n" if $strand eq "+";
      print join("\t", $scaffold, $genomic_end-1, $genomic_start), "\n" unless $strand eq "+";
      $start_exon_index++;
   }
}

close($bedfh);

exit 0;
