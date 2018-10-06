#!/usr/bin/perl
use POSIX;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $SCRIPTNAME = "subsetVCFstats.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

subsetVCFstats.pl - Subsets VCF stats from a TSV based on BED regions

=head1 SYNOPSIS

subsetVCFstats.pl [options]

 Options:
  --help,-h,-?         Display this help documentation
  --vcf_stats,-i       Input TSV of VCF statistics per site
                       e.g. derived from GATK VariantsToTable or
                       bcftools query
                       (default: -, which is STDIN)
  --bed_file,-b        BED file (3 columns) of regions to subset
  --debug,-d           Output extra information to STDERR

=head1 DESCRIPTION

subsetVCFstats.pl is designed to extract the lines of a TSV file
(assuming the first two columns are chromosome and position)
if they fall within the regions specified by the BED file.

This is a fairly generic operation, but this script was originally
written for processing variant calling simulation results.

=cut


my $stats_path = "";
my $bed_path = "";
my $help = 0;
my $man = 0;
my $debug = 0;
my $dispversion = 0;

GetOptions('vcf_stats|i=s' => \$stats_path, 'bed_file|b=s' => \$bed_path, 'debug|d+' => \$debug, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $bed_file;
open $bed_file, "<", $bed_path or die "Failed to open input BED file ${bed_path} due to error $!";

print STDERR "Reading BED file ${bed_path}\n" if $debug;
my %query_starts = ();
my %query_ends = ();
while (my $line = <$bed_file>) {
   chomp $line;
   my ($chrom, $bedstart, $end) = split /\t/, $line, 3;
   $query_starts{$chrom} = [] unless exists($query_starts{$chrom});
   $query_ends{$chrom} = [] unless exists($query_ends{$chrom});
   push @{$query_starts{$chrom}}, $bedstart+1;
   push @{$query_ends{$chrom}}, $end;
}
print STDERR "Found intervals for ", scalar(keys %query_starts), " scaffolds\n" if $debug;
close($bed_file);

print STDERR "Processing stats file ${stats_path}\n" if $debug;
my $stats_file;
if ($stats_path eq '-') {
   open $stats_file, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input stats TSV file due to error $!";
} else {
   open $stats_file, "<", $stats_path or die "Failed to open input VCF stats TSV file ${stats_path} due to error $!";
}
my $interval_index = 0;
my $prev_chrom = "";
while (my $line = <$stats_file>) {
   chomp $line;
   next if $line =~ /^#/; #Skip header line or any comment lines
   my ($chrom, $pos, $rest) = split /\t/, $line, 3;
   #Reset the interval index for each new scaffold:
   $interval_index = 0 unless $chrom eq $prev_chrom;
   $prev_chrom = $chrom;
   #Skip the site if there aren't any intervals for this scaffold:
   next unless exists($query_starts{$chrom});
   #Move through intervals until either contains or is just past pos:
   $interval_index++ until $interval_index == scalar(@{$query_starts{$chrom}}) or $query_ends{$chrom}[$interval_index] >= $pos;
   #Skip the site if there are no intervals left on this scaffold:
   next if $interval_index == scalar(@{$query_starts{$chrom}});
   #Skip the site if the interval starts after the current pos:
   next if $query_starts{$chrom}[$interval_index] > $pos;
   #Only print the site if we're within an interval:
   print $line, "\n" if $query_starts{$chrom}[$interval_index] <= $pos and $query_ends{$chrom}[$interval_index] >= $pos;
}
print STDERR "Done subsetting stats file ${stats_path}\n" if $debug;
close($stats_file);
