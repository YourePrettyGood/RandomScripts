#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.1 (2018/11/15) Abstracted revcomp to function      #
################################################################

#First pass script to extract CDSes identified in a file resembling BED,
# but with a fourth column identifying the protein to be extracted.

my $SCRIPTNAME = "extractCDSfromVulgar.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

extractCDSfromVulgar.pl - Extract CDSes from a genome given VULGAR alignments

=head1 SYNOPSIS

extractCDSfromVulgar.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --input_genome,-i     Path to input genome FASTA file (default: STDIN)
  --vulgar_path,-v      Path to filtered VULGAR alignment outputs from
                        exonerate (exclude the Exonerate head and tail)
  --header_prefix,-p    String to prepend to protein ID for output FASTA
                        (default is empty string, so _ is prepended)
  --debug,-d            Output debugging information to STDERR
  --version             Output version string

=head1 DESCRIPTION

This script extracts the CDSes indicated by a VULGAR alignment produced
by exonerate from the input FASTA file, and concatenates them into a
series of FASTA records on STDOUT.

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $vulgar_path = "";
my $header_prefix = "";
my $debug = 0;
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'vulgar_path|v=s' => \$vulgar_path, 'header_prefix|p=s' => \$header_prefix, 'version' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
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

#Open the exonerate alignment VULGAR file:
unless(open(VULGAR, "<", $vulgar_path)) {
   print STDERR "Error opening exonerate alignment VULGAR file.\n";
   exit 3;
}

sub vulgarToIntervals($$$) {
   my @vulgaroptuples = @{shift @_};
   my @intervals = ();
   my $num_optuples = scalar@vulgaroptuples/3;
   my $interval_start = shift @_;
   my $strand = shift @_;
   my $interval_end = $interval_start;
   my $previous_op = "M"; #Keep track of the previous operation, in case of adjacent skip operations
   my $shift_len; #Offset due to frameshift or gap operations
   for (my $i = 0; $i < $num_optuples; $i++) {
      my $op = $vulgaroptuples[3*$i];
      my $prot_len = $vulgaroptuples[3*$i+1];
      my $dna_len = $vulgaroptuples[3*$i+2];
      if ($op eq "M") {
         #$interval_end += $dna_len-1 if $strand eq "+";
         #$interval_end -= $dna_len-1 if $strand eq "-";
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
      } elsif ($op eq "S") {
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
      } elsif ($op eq "5" or $op eq "G" or $op eq "F") {
         #Only output the prior interval if you just were in a match region:
         #push @intervals, join("-", $interval_start, $interval_end) if $previous_op eq "M" or $previous_op eq "S";
         push @intervals, join("-", $interval_start, $interval_end-1) if ($previous_op eq "M" or $previous_op eq "S") and $strand eq "+";
         push @intervals, join("-", $interval_start, $interval_end+1) if ($previous_op eq "M" or $previous_op eq "S") and $strand eq "-";
         #$shift_len = $dna_len+1 if $previous_op eq "M" or $previous_op eq "S";
         #The only cases where the above is not true are adjacent F and/or G operations, where we only do the +1 once:
         #$shift_len = $dna_len unless $previous_op eq "M" or $previous_op eq "S";
         $shift_len = $dna_len; #Universal, now that we only shift the end back when outputting intervals
         $interval_end += $shift_len if $strand eq "+";
         $interval_end -= $shift_len if $strand eq "-";
         $interval_start = $interval_end;
      } elsif ($op eq "I" or $op eq "3") {
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
         $interval_start = $interval_end;
      } else {
         die "Unknown VULGAR op ${op}\n";
      }
      $previous_op = $op;
   }
   push @intervals, join("-", $interval_start, $interval_end-1) if $vulgaroptuples[-3] eq "M" and $strand eq "+";
   push @intervals, join("-", $interval_start, $interval_end+1) if $vulgaroptuples[-3] eq "M" and $strand eq "-";
   return(\@intervals);
}

my %CDS_intervals = ();
my %CDS_scaffolds = ();
my %CDS_strands = ();
while (my $line = <VULGAR>) {
   chomp $line;
   next unless $line =~ /^vulgar:/; #Make sure to skip non-VULGAR lines
   my @vulgararr = split /\s+/, $line;
   shift @vulgararr; #"vulgar:"
   my $prot_id = shift @vulgararr;
   shift @vulgararr; #Residue start
   shift @vulgararr; #Residue end
   shift @vulgararr; #Meaningless protein strand
   my $prot_scaf = shift @vulgararr;
   my $prot_start = shift @vulgararr;
   my $prot_end = shift @vulgararr;
   my $prot_strand = shift @vulgararr;
   $prot_start += 1 if $prot_strand eq "+";
   shift @vulgararr; #Alignment score
   $CDS_scaffolds{$prot_id} = $prot_scaf;
   $CDS_intervals{$prot_id} = vulgarToIntervals(\@vulgararr, $prot_start, $prot_strand);
   $CDS_strands{$prot_id} = $prot_strand;
}

close(VULGAR);

#Now we can iterate through the genome FASTA, and store the records
my %scaffolds = ();
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <GENOME>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the sites from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      $scaffolds{$scaffold_name} = $scaffold_sequence unless $scaffold_sequence eq "" or $scaffold_name eq "";
      my $scaffold_name_line = substr $line, 1; #Get rid of the prefixed ">"
      my @scaffold_name_parts = split /\s+/, $scaffold_name_line;
      $scaffold_name = $scaffold_name_parts[0];
      $scaffold_sequence = ""; #Clear out the old sequence
   } else { #Sequence line
      $scaffold_sequence .= $line;
   }
}
$scaffolds{$scaffold_name} = $scaffold_sequence;
close(GENOME);

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/AaCcGgTtRrYySsWwKkMmBbDdHhVvNn/TtGgCcAaYyRrSsWwMmKkVvHhDdBbNn/; #Complement incl. IUPAC degenerate bases
   return $reverse_sequence;
}

sub printCDS($$$$$) {
   my $prot_id = shift @_;
   my $header_prefix = shift @_;
   my $sequence = shift @_;
   my @intervals = @{shift @_;};
   my $strand = shift @_;
   my $CDS = "";
   for my $interval (@intervals) {
      my ($start, $end) = split /-/, $interval, 2;
      if ($strand eq "+") { # + strand
         my $segment = substr $sequence, $start-1, $end-$start+1;
         $CDS .= uc($segment);
      } else { # - strand
         my $segment = substr $sequence, $end-1, $start-$end+1;
         $CDS .= uc(revcomp($segment));
      }
   }
   print STDOUT ">", $header_prefix, "_", $prot_id, "\n", $CDS, "\n";
}

for my $CDS (keys %CDS_scaffolds) {
   my $scaffold = $CDS_scaffolds{$CDS};
   my $intervals = $CDS_intervals{$CDS};
   my $strand = $CDS_strands{$CDS};
   die "Scaffold ${scaffold} not found in input FASTA.\n" unless exists($scaffolds{$scaffold});
   my $scaf_seq = $scaffolds{$scaffold};
   print STDERR $CDS, " ", $scaffold, " ", join(":", @{$intervals}), "\n" if $debug; #Diagnostic
   printCDS($CDS, $header_prefix, $scaf_seq, $intervals, $strand);
}

exit 0;
