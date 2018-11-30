#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

###############################################################################
#                                                                             #
# Version 1.1 (2018/11/22) Usable fraction for col. 4 (compatibility)         #
###############################################################################

#Add records for sites not output by other programs
#For example, sites not covered by 1:1 LAST alignments are
#not output by divergenceFromMAF.pl

my $SCRIPTNAME = "decompressStats.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

decompressStats.pl - Add missing records to genome-wide per-site stats output

=head1 SYNOPSIS

decompressStats.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --input_stats,-i      Path to input stats file (default: STDIN)
  --fai,-f              Path to FASTA index for genome assembly
  --stat_column,-s      Which column contains the statistic? (default: 3)
  --usable_fraction,-u  Output usable fraction instead of omit as col. 4
  --debug,-d            Output debugging information to STDERR
  --version,-v          Display the version string

=head1 DESCRIPTION

This script adds missing records to genome-wide per-site stats
files output by other programs (following the 4-column format
Scaffold ID, Position, Stat, Omit). Deviations from the 4-column
format are permitted if you specify which column in the input
contains the statistic of interest. Missing records are assumed
to be requiring omission. This script is generally piped in
between the other program and a windowing/summary program
like nonOverlappingWindows.

For example, given a 2L chromosome of length 15 and stats output as follows:

=begin text

2L	10	1	0
2L	11	0	0
2L	12	0	0
2L	13	1	0
2L	14	0	0
2L	15	0	0

=end text

this script should output:

=begin text

2L	1	0	1
2L	2	0	1
2L	3	0	1
2L	4	0	1
2L	5	0	1
2L	6	0	1
2L	7	0	1
2L	8	0	1
2L	9	0	1
2L	10	1	0
2L	11	0	0
2L	12	0	0
2L	13	1	0
2L	14	0	0
2L	15	0	0

=end text

=cut

my $help = 0;
my $man = 0;
my $stats_path = "STDIN";
my $fai_path = "";
my $stat_column = 3;
my $usable_fraction = 0;
my $debug = 0;
my $dispversion = 0;
GetOptions('input_stats|i=s' => \$stats_path, 'fai|f=s' => \$fai_path, 'stat_column|s=i' => \$stat_column, 'usable_fraction|u' => \$usable_fraction, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Check that the stat column is not 1, 2, or 4:
if ($stat_column == 1 or $stat_column == 2 or $stat_column == 4) {
   print STDERR "Invalid choice for statistic column: ${stat_column}\n";
   exit 5;
}

#Open the stats file, or set it up to be read from STDIN:
print STDERR "Opening stats file\n" if $debug;
my $stats_fh;
if ($stats_path ne "STDIN") {
   unless(open($stats_fh, "<", $stats_path)) {
      print STDERR "Error opening stats file: ${stats_path}.\n";
      exit 2;
   }
} else {
   open($stats_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $stats_fh so we can seamlessly handle piping
}

my $fai_fh;
if (-e $fai_path) {
   unless(open($fai_fh, "<", $fai_path)) {
      print STDERR "Error opening fai file: ${fai_path}.\n";
      exit 3;
   }
} else {
   print STDERR "Provided fai file ${fai_path} does not exist.\n";
   exit 4;
}

#Read in the fai so that we know the scaffold boundaries:
print STDERR "Reading in fai file ${fai_path}\n" if $debug;
my %scaffold_length = ();
my %scaffold_offset = ();
while (my $line = <$fai_fh>) {
   chomp $line;
   my @fai_elems = split /\t/, $line;
   $scaffold_length{$fai_elems[0]} = $fai_elems[1];
   $scaffold_offset{$fai_elems[0]} = 1;
}
close($fai_fh);

print STDERR "Done processing ", scalar(keys %scaffold_length), " scaffolds from fai file ${fai_path}\n" if $debug;

print STDERR "Processing stats file ${stats_path}\n" if $debug;
my $prev_scaffold = "";
my %used_scaffolds = ();
my $empty_extra_columns = "";
while (my $line = <$stats_fh>) {
   chomp $line;
   my @stats_line = split /\t/, $line;
   my $scaffold = $stats_line[0];
   my $position = $stats_line[1];
   my $omit_position = $stats_line[3];
   my $num_columns = scalar(@stats_line);
   if ($stat_column > $num_columns) {
      print STDERR "Selected statistic column (${stat_column}) does not exist in the input ${stats_path}\n";
      close($stats_fh);
      exit 6;
   }
   my $stat_value = $stats_line[$stat_column-1];
   $empty_extra_columns = "\t0"x(${num_columns}-4);
   $used_scaffolds{$scaffold} = 1;
   if ($scaffold ne $prev_scaffold and $prev_scaffold ne "") {
      print STDERR "Done processing scaffold ${prev_scaffold}\n" if $debug;
      while ($scaffold_offset{$prev_scaffold} <= $scaffold_length{$prev_scaffold}) {
         print STDOUT join("\t", $prev_scaffold, $scaffold_offset{$prev_scaffold}, 0, 1-$usable_fraction), $empty_extra_columns, "\n";
         $scaffold_offset{$prev_scaffold}++;
      }
   }
   while ($scaffold_offset{$scaffold} < $position) {
      print STDOUT join("\t", $scaffold, $scaffold_offset{$scaffold}, 0, 1-$usable_fraction), $empty_extra_columns, "\n";
      $scaffold_offset{$scaffold}++;
   }
   print STDOUT $line, "\n";
   $scaffold_offset{$scaffold}++;
   $prev_scaffold = $scaffold;
}
print STDERR "Done processing stats file ${stats_path}\n" if $debug;
close($stats_fh);

#Fill out the last scaffold if necessary:
print STDERR "Filling in missing records for the last scaffold\n" if $debug;
while ($scaffold_offset{$prev_scaffold} <= $scaffold_length{$prev_scaffold}) {
   print STDOUT join("\t", $prev_scaffold, $scaffold_offset{$prev_scaffold}, 0, 1-$usable_fraction), $empty_extra_columns, "\n";
   $scaffold_offset{$prev_scaffold}++;
}

#Fill out missing scaffolds if necessary:
print STDERR "Filling in missing records for omitted scaffolds\n" if $debug;
for my $scaffold (keys %scaffold_length) {
   unless (exists($used_scaffolds{$scaffold})) {
      while ($scaffold_offset{$scaffold} <= $scaffold_length{$scaffold}) {
         print STDOUT join("\t", $scaffold, $scaffold_offset{$scaffold}, 0, 1-$usable_fraction), $empty_extra_columns, "\n";
         $scaffold_offset{$scaffold}++;
      }
   }
}
print STDERR "Process complete!" if $debug;

exit 0;
