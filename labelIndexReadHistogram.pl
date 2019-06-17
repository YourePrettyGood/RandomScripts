#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.1 (2019/04/12) Extra bit of logging                #
################################################################

#First-pass script to label barcodes in an index read histogram
# based on a barcode file, and incorporating mismatches if no
# perfectly matching barcode is found.

my $SCRIPTNAME = "labelIndexReadHistogram.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

labelIndexReadHistogram.pl - Label an index read histogram with sample IDs

=head1 SYNOPSIS

labelIndexReadHistogram.pl [options]

 Options:
  --help,-h,-?         Display this documentation
  --input_histogram,-i Path to input index read histogram TSV file (default: STDIN)
  --barcode_file,-b    Path to barcode file
  --version,-v         Output version string

=head1 DESCRIPTION

Using the information found in the barcode file, this script adds a third
column to the index read histogram that indicates the sample ID and how
many mismatches were required for the barcode sequence to match.

=cut

my $help = 0;
my $man = 0;
my $histogram_path = "STDIN";
my $bcfile_path = "";
my $dispversion = 0;
GetOptions('input_histogram|i=s' => \$histogram_path, 'barcode_file|b=s' => \$bcfile_path, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the index read histogram file, or set it up to be read from STDIN:
my $histfile_fh;
if ($histogram_path ne "STDIN") {
   unless(open($histfile_fh, "<", $histogram_path)) {
      print STDERR "Error opening index read histogram TSV file ${histogram_path}.\n";
      exit 1;
   }
} else {
   open($histfile_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $histfile_fh so we can seamlessly handle piping
}

#Open the barcode file:
my $indexfile_fh;
unless(open($indexfile_fh, "<", $bcfile_path)) {
   print STDERR "Error opening barcode file ${bcfile_path}.\n";
   exit 2;
}

my @bases = ("A", "C", "G", "T");
my %FCmap = ();
while (my $line = <$indexfile_fh>) {
   chomp $line;
   my ($id, $bc) = split /\t/, $line, 2;
   #Should do some case checking to make sure it's a valid barcode file in the future
   die "No tabs present in barcode file, make sure your columns are separated by tabs!\n" unless defined($bc);
   die "Too many columns in barcode file, expecting 2.\n" if $bc =~ /\t/;
   die "Invalid barcode found, must be non-ambiguous nucleotides (ACGT only).\n" unless $bc =~ /^[ACGT]+$/;
   $FCmap{$bc} = $id;
}
close($indexfile_fh);

while (my $line = <$histfile_fh>) {
   chomp $line;
   $line =~ s/^\s+//;
   my ($count, $bc) = split /\s+/, $line, 2;
   print $count, "\t", $bc, "\t";
   if (exists($FCmap{$bc})) {
      print $FCmap{$bc}, " (0 mismatches)";
   } else {
      for (my $pos = 0; $pos < length($bc); $pos++) {
         for my $base (@bases) {
            my $altbc = $bc;
            substr($altbc, $pos, 1) = $base;
            next if $altbc eq $bc;
            if (exists($FCmap{$altbc})) {
               print $FCmap{$altbc}, " (1 mismatch)";
               last;
            }
            for (my $postwo = $pos+1; $postwo < length($bc); $postwo++) {
               for my $basetwo (@bases) {
                  my $twoaltbc = $altbc;
                  substr($twoaltbc, $postwo, 1) = $basetwo;
                  next if $twoaltbc eq $bc or $twoaltbc eq $altbc;
                  if (exists($FCmap{$twoaltbc})) {
                     print $FCmap{$twoaltbc}, " (2 mismatches)\t";
                  }
               }
            }
         }
      }
   }
   print "\n";
}
close($histfile_fh);
