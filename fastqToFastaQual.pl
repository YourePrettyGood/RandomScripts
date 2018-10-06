#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#########################################################################################
# fastqToFastaQual.pl                                                                   #
# Version 1.0 (2017/04/20)                                                              #
# Description:                                                                          #
# This script converts a FASTQ file into a paired set of FASTA and QUAL files           #
#                                                                                       #
# Usage:                                                                                #
#  fastqToFastaQual.pl [-i input.fastq] [-a output.fasta] [-q output.qual]              #
# Options:                                                                              #
#  --input_file,-i:   		Input FASTQ file name (default: STDIN)                  #
#  --fasta_file,-a:   		Output FASTA file name                                  #
#  --qual_file,-q:		Output QUAL file name                                   #
#########################################################################################

my $SCRIPTNAME = "fastqToFastaQual.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

fastqToFastaQual.pl - Generate FASTA and QUAL from FASTQ for GenBank submission

=head1 SYNOPSIS

fastqToFastaQual.pl [options]

 Options:
  --help,-h,-?		Display this help documentation
  --input_file,-i	Input FASTQ file name (default: STDIN)
  --fasta_file,-a	Output FASTA file
  --qual_file,-q	Output QUAL file
  --version,-v          Output version string

=head1 DESCRIPTION
This script generates FASTA and QUAL files from a FASTQ file, which may be
useful for Genbank submissions from Quivered assemblies.

=cut

my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $output_fasta = "";
my $output_qual = "";
my $dispversion = 0;
GetOptions('input_file|i=s' => \$input_path, 'fasta_file|a=s' => \$output_fasta, 'qual_file|q=s' => \$output_qual, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

if ($input_path ne "STDIN") {
   unless(open(FASTQ, "<", $input_path)) {
      print STDERR "Error opening input FASTQ file.\n";
      exit 2;
   }
} else {
   open(FASTQ, "<&", "STDIN"); #Duplicate the file handle for STDIN to CONTIGS so we can seamlessly handle piping
}

if ($output_fasta eq "") {
   print STDERR "Missing FASTA file path for output.\n";
   exit 5;
} else {
   unless(open(FASTA, ">", $output_fasta)) {
      print STDERR "Error opening output FASTA file.\n";
      exit 4;
   }
}

if ($output_qual eq "") {
   print STDERR "Missing QUAL file path for output.\n";
   exit 7;
} else {
   unless(open(QUAL, ">", $output_qual)) {
      print STDERR "Error opening output QUAL file.\n";
      exit 6;
   }
}

sub wrap_fasta_sequence($$) {
   my $fasta_sequence = shift @_;
   my $line_length = shift @_;
   my $sequence_length = length($fasta_sequence);
   my @wrapped_lines = ();
   for (my $line_start = 0; $line_start < $sequence_length; $line_start += $line_length) {
      my $wrapped_line = substr($fasta_sequence, $line_start, $line_length);
      $wrapped_line =~ s/^\s+//; #Remove leading whitespace
      $wrapped_line =~ s/\s+$//; #Remove trailing whitespace
      push @wrapped_lines, $wrapped_line;
   }
   return(join("\n", @wrapped_lines) . "\n");
}

sub wrap_qual_sequence($$) {
   my $qual_sequence = shift @_;
   my $line_length = shift @_;
   #die "Invalid qual line wrapping length -- must give remainder 2 when dividing by 3\n" unless $line_length % 3 == 2;
   my @quals = ();
   my @qualarr = split //, $qual_sequence;
   for my $qual (@qualarr) {
      push @quals, sprintf("%02d", ord($qual)-33);
   }
   return(wrap_fasta_sequence(join(" ", @quals), $line_length));
}

my ($header, $sequence, $qualheader, $quals) = ('', '', '', '');
my $fastq_line_modulus = 1;
while (my $line = <FASTQ>) {
   chomp $line;
   if ($fastq_line_modulus == 1) {
      $header = substr $line, 1;
   } elsif ($fastq_line_modulus == 2) {
      $sequence = $line;
      my $wrapped_sequence = wrap_fasta_sequence($sequence, 80);
      print FASTA ">", $header, "\n", $wrapped_sequence, "\n";
   } elsif ($fastq_line_modulus == 3) {
      $qualheader = substr $line, 1;
   } else { #Modulus is 0
      $quals = $line;
      my $wrapped_quals = wrap_qual_sequence($quals, 51);
      print QUAL ">", $header, "\n", $wrapped_quals, "\n";
      ($header, $sequence, $qualheader, $quals) = ('', '', '', '');
   }
   $fastq_line_modulus++;
   $fastq_line_modulus %= 4;
}
#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close(FASTQ);
}
close(FASTA);
close(QUAL);

exit 0;
