#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.1 (2018/11/15) Revcomp now works on IUPAC bases    #
# Version 1.2 (2019/04/11) Non-greedy Parent regex and prefix  #
# Version 1.3 (2019/04/27) Option for longest isoform only     #
# Version 1.4 (2019/05/16) Stable choice of longest isoform    #
# Version 1.5 (2019/06/13) Use exons instead of CDS on request #
################################################################

#First pass script to construct a FASTA of CDSes from a GFF3 and
# a genome FASTA.  This uses principles from fixExons.pl for
# constructing the full CDS from an Exon Range String, and
# computes the Exon Range Strings from the GFF3.

#Revised 2018/01/15 to deal with GFF3 containing more than just
# CDS records, although still assuming local sortedness
# (i.e. sorted order within a gene)

#2019/05/16 revision just forces a sort order for the longest
# isoform choice so that the same transcript is chosen if there
# is a tie for longest.

#2019/06/13 revision allows for defining transcripts with
# child tags other than CDS (e.g. exon)

my $SCRIPTNAME = "constructCDSesFromGFF3.pl";
my $VERSION = "1.5";

=pod

=head1 NAME

constructCDSesFromGFF3.pl - Construct a FASTA of CDSes based on a GFF3 and genome FASTA

=head1 SYNOPSIS

constructCDSesFromGFF3.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_genome,-i      Path to input genome FASTA file (default: STDIN)
  --gff3_file,-g         Path to genome annotation GFF3 file
  --prefix,-p            Prefix to prepend to FASTA headers (default: none)
  --output_ers,-e        Output the "Exon Range String" for the CDS?
                         0 or 1, default is 0
  --longest,-l           Only output the longest isoform for each gene
  --feature,-f           Name of GFF3 feature to splice (default: CDS)
  --version,-v           Output version string

=head1 DESCRIPTION

This script constructs CDSes based on an input genome FASTA file and a
corresponding genome annotation GFF3 file.  The output is in a FASTA-like
format (no wrapping).

Optionally, it can be used to construct full transcripts (including UTRs)
as long as a single feature name represents all the components of a full
transcript.

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
my $header_prefix = "";
my $output_exon_range_strings = 0;
my $longest_only = 0;
my $feature_type = "CDS";
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'gff3_file|g=s' => \$gff3_path, 'prefix|p=s' => \$header_prefix, 'output_ers|e' => \$output_exon_range_strings, 'longest|l' => \$longest_only, 'feature|f=s' => \$feature_type, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;
$header_prefix .= "_" unless $header_prefix eq "";

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

print STDERR "Unsupported feature type ${feature_type}, please use CDS or exon\n" unless $feature_type eq "CDS" or $feature_type eq "exon";
exit 4 unless $feature_type eq "CDS" or $feature_type eq "exon";

#Open the genome FASTA file, or set it up to be read from STDIN:
my $genome_fh;
if ($genome_path ne "STDIN") {
   unless(open($genome_fh, "<", $genome_path)) {
      print STDERR "Error opening genome FASTA file ${genome_path}.\n";
      exit 2;
   }
} else {
   open($genome_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to GENOME so we can seamlessly handle piping
}

#Open the annotation GFF3 file, which should solely contain sorted CDS records:
my $gff_fh;
unless(open($gff_fh, "<", $gff3_path)) {
   print STDERR "Error opening GFF3 file ${gff3_path}.\n";
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

my $FASTA_skip = 0; #Skip all lines after the ##FASTA line if present

#Iterate over the GFF features and build the Exon Range Strings
# for each CDS:
my %CDS_coordinates = ();
my %CDSes_per_scaffold = ();
my %isoform_length = (); #Keys are transcript IDs, values are isoform lengths
my %gene_tx_map = (); #Hash of hashes ({gene}{tx}), basically a set of sets
while (my $line = <$gff_fh>) {
   chomp $line;
   $FASTA_skip = 1 if $line =~ /^##FASTA/; #Trigger FASTA skip
   next if $FASTA_skip;
   next if $line =~ /^#/; #Skip comment lines
   my ($scaffold, $set, $type, $start, $end, $score, $strand, $frame, $tag_string) = split /\t/, $line, 9;
   if ($type =~ /mRNA/i) {
      if ($tag_string =~ /Parent=(.+?)(?:;|$)/i) {
         my $gene_ID = $1;
         if ($tag_string =~ /ID=(.+?)(?:;|$)/i) {
            my $tx_ID = $1;
            $gene_tx_map{$gene_ID} = {} unless exists($gene_tx_map{$gene_ID});
            $gene_tx_map{$gene_ID}{$tx_ID} = undef;
            $isoform_length{$tx_ID} = 0; #Initialize the isoform length sum
         } else {
            print STDERR "Regex to find transcript ID failed for tag string: ", $tag_string, "\n";
            next;
         }
      } else {
         print STDERR "Regex to find gene ID failed for tag string: ", $tag_string, "\n";
         next;
      }
   }
   next unless $type eq $feature_type;
   my @transcript_names = ();
   if ($tag_string =~ /Parent=(.+?)(?:;|$)/i) {
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
      #Add this CDS' length to the isoform length sum:
      $isoform_length{$transcript_name} += $end - $start + 1;
   }
}

close($gff_fh);

#Determine the longest isoform for each gene:
my %longest_isoform = (); #Keys are transcript IDs, basically a set
for my $gene (sort keys %gene_tx_map) {
   my $max_len = 0;
   my $max_len_id = "";
   for my $tx (sort keys %{$gene_tx_map{$gene}}) {
      if ($isoform_length{$tx} > $max_len) {
         $max_len = $isoform_length{$tx};
         $max_len_id = $tx;
      }
   }
   $longest_isoform{$max_len_id} = undef;
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
while (my $line = <$genome_fh>) {
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
            next if $longest_only and !exists($longest_isoform{$CDS});
            #Output the FASTA record of the CDS to STDOUT:
            print STDOUT ">", $header_prefix, $CDS, "\n", $constructed_CDS, "\n";
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
close($genome_fh);
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
      unless ($longest_only and !exists($longest_isoform{$CDS})) {
         #Output the FASTA record of the CDS to STDOUT:
         print STDOUT ">", $header_prefix, $CDS, "\n", $constructed_CDS, "\n";
         print STDOUT $scaffold_name, "=", $exon_range_string, "\n" if $output_exon_range_strings;
         print STDOUT $scaffold_name, "=", computeCDSRangeString($exon_range_string), "\n" if $output_exon_range_strings;
      }
   }
}

exit 0;
