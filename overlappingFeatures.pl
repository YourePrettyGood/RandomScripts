#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

my $SCRIPTNAME = "overlappingFeatures.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

overlappingFeatures.pl - Identify GFF features overlapping BED intervals

=head1 SYNOPSIS

overlappingFeatures.pl [options]

 Options:
  --help,-h,-?          Display this help documentation
  --interval_bed,-b     Input BED file indicating intervals of interest
                        (Default: STDIN)
  --input_gff,-g        Input GFF3 file (sorted) of features for the genome
  --version,-v          Output version string

=head1 DESCRIPTION

overlappingFeatures.pl outputs the GFF3 lines (plus an extra column of the
BED interval of interest overlapping the GFF3 feature(s)) that overlap with
intervals of interest provided in BED 3-column format.  The idea is fairly
simple, but essentially allows identification of features underlying
potentially interesting base-level signals (e.g. common IBD tracts, high
Dxy, etc.).

=cut

#Parse options:
my $help = 0;
my $man = 0;
my $bed_path = '';
my $gff_path = '';
my $dispversion = 0;
GetOptions('input_gff|g=s' => \$gff_path, 'interval_bed|b=s' => \$bed_path, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my %intervals = ();

my $bed;
if ($bed_path eq '' or $bed_path eq '-') {
   open $bed, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input BED due to error $!\n";
} else {
   open $bed, "<", $bed_path or die "Failed to open input BED file due to error $!\n";
}
while (my $line = <$bed>) {
   chomp $line;
   my ($scaf, $bedstart, $bedend) = split /\t/, $line, 3;
   $intervals{$scaf} = () unless exists($intervals{$scaf});
   push @{$intervals{$scaf}}, join("-", $bedstart+1, $bedend);
}
close($bed);

my $gff;
open($gff, "<", $gff_path) or die "Unable to open GFF file $gff_path due to error $!";
while (my $line = <$gff>) {
   chomp $line;
   my @gfflinearr = split /\t/, $line;
   for my $interval (@{$intervals{$gfflinearr[0]}}) {
      my ($intstart, $intend) = split /-/, $interval, 2;
      print $line, "\t", $interval, "\n" if $gfflinearr[3] < $intend and $gfflinearr[4] > $intstart;
      last if $gfflinearr[3] < $intend and $gfflinearr[4] > $intstart;
   }
}
close($gff);
