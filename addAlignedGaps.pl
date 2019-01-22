#!/usr/bin/env perl

use warnings;
use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

######################################################################
# addAlignedGaps.pl                                                  #
# Usage:                                                             #
#  addAlignedGaps.pl [-a input alignment] [-i unaligned FASTAs]      #
#    [-o output aligned FASTA]                                       #
#                                                                    #
# Arguments:                                                         #
#  -a,--input_alignment   Aligned multi-FASTA of references where    #
#                         the prefix up to _ is the species name     #
#                         (mandatory, no default)                    #
#  -i,--input_FASTA       FASTA of pseudoreferences to put into the  #
#                         aligned coordinate space, where prefixes   #
#                         up to _ are the species name, matching one #
#                         of the species in the alignment            #
#                         (default: STDIN)                           #
#  -o,--output_alignment  Filename for output aligned multi-FASTA    #
#                         containing sequences from both the input   #
#                         alignment and gap-inserted versions of the #
#                         unaligned sequences                        #
#                         (default: STDOUT)                          #
# Description:                                                       #
#  addAlignedGaps.pl produces an aligned FASTA of all sequences based#
#  on adding gaps found in an alignment of reference sequences to a  #
#  set of unaligned sequences, each set in an identifiable reference #
#  coordinate space. We use the prefix of the FASTA header before the#
#  first _ to identify the reference, and to match to equivalent     #
#  prefixes in the unaligned FASTA headers.                          #
######################################################################

my $SCRIPTNAME = "addAlignedGaps.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

addAlignedGaps.pl - Add gaps to unaligned sequences and append to MSA

=head1 SYNOPSIS

addAlignedGaps.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  -a,--input_alignment   Aligned multi-FASTA of references where
                         the prefix up to _ is the species name
                         (mandatory, no default)
  -i,--input_FASTA       FASTA of pseudoreferences to put into the
                         aligned coordinate space, where prefixes
                         up to _ are the species name, matching one
                         of the species in the alignment
                         (default: STDIN)
  -o,--output_alignment  Filename for output aligned multi-FASTA
                         containing sequences from both the input
                         alignment and gap-inserted versions of the
                         unaligned sequences
                         (default: STDOUT)
  --version,-v           Output version string
  --debug,-d             Print diagnostic output to STDERR

=head1 DESCRIPTION

addAlignedGaps.pl produces an aligned FASTA of all sequences based
on adding gaps found in an alignment of reference sequences to a
set of unaligned sequences, each set in an identifiable reference
coordinate space. We use the prefix of the FASTA header before the
first _ to identify the reference, and to match to equivalent
prefixes in the unaligned FASTA headers.

=cut

sub gapPositions($) {
   my $seq = shift @_;
   my @gap_positions = ();
   my @bases = split //, $seq;
   my $seqlen = scalar(@bases);
   for (my $i = 0; $i < $seqlen; $i++) {
      if ($bases[$i] eq "-") {
         push @gap_positions, $i;
      }
   }
   return \@gap_positions;
}

sub addGaps($$$) {
   my $seq = shift @_;
   my $species = shift @_;
   my %gap_sequences = %{shift @_;};
   my @gap_sequence = ();
   if (exists($gap_sequences{$species})) {
      @gap_sequence = @{$gap_sequences{$species}};
   } else {
      print STDERR "Species ${species} found in unaligned FASTA but not in alignment\n";
      return undef;
   }
   my @gapped_sequence = split //, $seq; #We're going to splice in gaps backwards
   my $gap_index = $#gap_sequence;
   while ($gap_index >= 0) {
      splice(@gapped_sequence, $gapped_sequence[$gap_index], 0, '-');
      #TODO: Might be somewhat faster to RLE the gap sequences
      $gap_index--;
   }
   return join("", @gapped_sequence);
}

#Parse options:
my $help = 0;
my $man = 0;
my $aln_path = '';
my $in_path = '';
my $out_path = '';
my $dispversion = 0;
my $debug = 0;
GetOptions('input_alignment|a=s' => \$aln_path, 'input_FASTA|i=s' => \$in_path, 'output_FASTA|o=s' => \$out_path, 'version|v' => \$dispversion, 'help|h|?+' => \$help, 'debug|d+' => \$debug, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open input alignment and output files:
my $aln;
if (! -e $aln_path) {
   print STDERR "Unable to find input aligned FASTA ${aln_path}, exiting\n";
   exit 2;
}
open $aln, "<", $aln_path or die "Failed to open aligned FASTA file ${aln_path} due to error $!\n";
my $out;
if ($out_path =~ /\.gz$/) {
   $out = new IO::Compress::Gzip $out_path or die "Failed to create Gzipped output FASTA file $out_path due to error $GzipError\n";
} elsif ($out_path eq '') {
   open $out, ">&", \*STDOUT or die "Failed to duplicate STDOUT file handle for output FASTA due to error $!\n";
} else {
   open $out, ">", $out_path or die "Failed to open output FASTA file ${out_path} due to error $!\n";
}

#Read in the necessary gap sequences from the alignment:
print STDERR "Reading in gap sequences from alignment\n" if $debug;
my %gap_sequences = ();
my $species = "";
my $gap_sequence = "";
while (my $line = <$aln>) {
   print $out $line; #Feed through to the output before processing
   if ($line =~ /^>/) { #Header line
      if ($gap_sequence ne "" and $species ne "") {
         $gap_sequences{$species} = gapPositions($gap_sequence);
         $gap_sequence = "";
      }
      if ($line =~ /^>([A-Za-z0-9]+)_/) {
         $species = $1; #Store the species name/ID via regex capture
      } else {
         print STDERR "Unable to detect an appropriate species name in FASTA header ${line}\nExpecting alphanumeric sequence before _\n";
         close $aln;
         close $out;
         exit 3;
      }
   } else {
      chomp $line;
      $gap_sequence .= $line; #Add sequence to buffer
   }
}
if ($gap_sequence ne "" and $species ne "") {
   @{$gap_sequences{$species}} = gapPositions($gap_sequence);
   $gap_sequence = "";
}
close $aln;
print STDERR "Done reading in gap sequences from alignment\n" if $debug;

#Open up the unaligned sequence FASTA for streaming:
print STDERR "Adding gaps to unaligned sequences\n" if $debug;
my $in;
if ($in_path =~ /\.gz$/) {
   $in = new IO::Uncompress::Gunzip $in_path or die "Failed to open Gzipped input FASTA file $in_path due to error $GunzipError\n";
} elsif ($in_path eq '') {
   open $in, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input FASTA due to error $!\n";
} else {
   open $in, "<", $in_path or die "Failed to open input FASTA file ${in_path} due to error $!\n";
}

#Iterate over the scaffolds:
$species = ""; #Keep track of the species prefix to extract the gap sequence
my $sequence = ""; #Store sequence progressively for wrapped FASTAs
while (my $line = <$in>) {
   if ($line =~ />/) {
      if ($sequence ne "" and $species ne "") {
         my $gapped_sequence = addGaps($sequence, $species, \%gap_sequences);
         unless (defined($gapped_sequence)) { #Catch the exception when species isn't in the hash
            close $in;
            close $out;
            exit 5;
         }
         print $out $gapped_sequence, "\n";
         $sequence = "";
      }
      print $out $line;
      if ($line =~ /^>([A-Za-z0-9]+)_/) {
         $species = $1; #Store the species name/ID via regex capture
      } else {
         print STDERR "Unable to detect an appropriate species name in FASTA header ${line}\nExpecting alphanumeric sequence before _\n";
         close $in;
         close $out;
         exit 4;
      }
   } else {
      chomp $line;
      $sequence .= $line; #Add line to the sequence buffer
   }
}
#Make sure to catch the final sequence:
if ($sequence ne "" and $species ne "") {
   my $gapped_sequence = addGaps($sequence, $species, \%gap_sequences);
   unless (defined($gapped_sequence)) { #Catch the exception when species isn't in the hash
      close $in;
      close $out;
      exit 5;
   }
   print $out $gapped_sequence, "\n";
   $sequence = "";
}


#Close input and output files:
close $in;
close $out;

exit 0;
