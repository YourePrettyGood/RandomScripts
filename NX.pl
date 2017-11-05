#!/usr/bin/perl -w
use POSIX;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

=pod

=head1 NAME

NX.pl - Calculates length-weighted quantile statistics from FASTAs

=head1 SYNOPSIS

NX.pl [options] <FASTA> [X] [G]

 Positional arguments:
  <FASTA>              Input FASTA file (may use - for STDIN)
  [X]                  Quantile to use (may specify as option instead)
                       (Default: 50)
  [G]                  Genome size to use instead of assembly size
                       (May specify as option instead)

 Options:
  --help,-h,-?         Display this help documentation
  --genome_size,-g     Genome size to use instead of assembly size
                       (Or specify as 3rd positional argument)
                       (Calculates an NGX instead of NX)
  --genome_percent,-p  Percentage of genome to use for quantile
                       (Or specify as 2nd positional argument)
                       (Default: 50)
  --debug,-d           Output extra information to STDERR

=head1 DESCRIPTION

NX.pl is a generalized script for calculating length-weighted quantile
statistics from a FASTA file.  The original script calculated N50,
but has since been adapted to calculate any NX or NGX statistic, like
N10, N90, NG50, etc.  It also outputs some baseline statistics like
the total assembly size, longest and shortest scaffold/contig, and
the average scaffold/contig size, plus the total number of scaffolds/
contigs in the assembly.

Hypothetically, you could also calculate read set statistics with NX.pl,
like the NR25 of a dataset by pre-calculating the number of bases for
25x coverage of the genome, specifying that as the genome size, and
using X=100.  This may take a while though, due to sorting.

=cut

my $input_filepath = '';
my $X = 50;
my $customsum = 0;
my $seqlen;
my $numseqs = 0;
my $help = 0;
my $man = 0;
my $debug = 0;

#This script calculates the NX of a set of contigs in FASTA format, where X is in (0,100].
#For X=50, this script is equivalent to N50.pl
#If a third argument is entered, this script uses that number instead of the total sum of contig lengths.
#Thus, the NG50 is calculated as ./NX.pl [FASTA file] 50 [genome size], and the NR25 is calculated as
#./NX.pl [read file] 100 [25*genome size]

GetOptions('genome_size|g=i' => \$customsum, 'genome_percent|p=i' => \$X, 'debug|d' => \$debug, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;


if ( scalar@ARGV > 0) { #Retrieve contig file path from command line argument list
   $input_filepath = $ARGV[0];
   if (scalar@ARGV > 1) { #Get the X (default 50)
      $X = $ARGV[1] + 0; #Make sure it's read as numeric
      die "Invalid X: Should be in the interval (0,100].\n" if $X <= 0 or $X > 100;
      if (scalar@ARGV > 2) { #Get the contig sum replacement, e.g. the G in NG50
         $customsum = $ARGV[2] + 0; #Make sure it's read as numeric
         die "Invalid custom threshold: Should be an integer >= 1.\n" if $customsum <= 0;
      }
   }
} else {
   print "Missing input FASTA file\n";
   exit;
}

my $contigfile;
if ($input_filepath eq '-') {
   open $contigfile, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input FASTA due to error $!\n";
} elsif ($input_filepath =~ /\.gz$/) {
   $contigfile = new IO::Uncompress::Gunzip $input_filepath or die "Failed to open input gzipped FASTA file due to error ${GunzipError}\n";
} else {
   open $contigfile, "<", $input_filepath or die "Failed to open input FASTA file due to error $!\n";
}
print STDERR "Now reading contig file\n" if $debug;
#Much of this code is based on subprogramQDD.pm by Dr. Emese Meglecz
my %contiglens = (); #Associative array (aka hash) to hold the contig length counts
my @contiglensarr; #Numerically-indexed array (aka list) to hold the contig lengths
my $contiglensum = 0;
$/ = ">"; #Change the line delimiter to > so that each FASTA record is read in at once
while (my $record = <$contigfile>) {
   $record =~ s/>//; #Remove the > prefix of the header line of the record
   unless ($record eq '') { #Make sure the record isn't empty
      $record =~ /.*\n/; #Match the pattern (i.e. everything before the newline) within $record's contents
      my $header = $&; #Copy the previous pattern match (i.e. everything before the newline) into $header
      $record =~ s/$header//; #Delete (i.e. substitute with nothing) the header from $record
      $record =~ s/\s//gm; #Delete all newlines and "whitespace characters" from the sequence
      #Now we have the sequence in $record, and the header in $header
      $seqlen = length $record;
      #Add the contig length to the running sum
      $contiglensum = $contiglensum + $seqlen;
      #Increment the number of contigs:
      $numseqs++;
      #Increment the number of contigs with that length
      if (exists $contiglens{$seqlen}) {
         $contiglens{$seqlen}++;
      } else {
         $contiglens{$seqlen} = 1;
         push(@contiglensarr, $seqlen);
      }
   }
}
close $contigfile;
print STDERR "Finished reading contig file\n" if $debug;
print "Total sum of lengths of contigs is $contiglensum\n";
print STDERR "Now sorting the contig length array\n" if $debug;
my $arrlen = scalar@contiglensarr;
print STDERR "Contig length array has $arrlen elements\n" if $debug;

my @sortedcontiglens = sort {$b <=> $a} @contiglensarr;
print STDERR "Finished sorting the contig length array\n" if $debug;
#Now progressively add contigs from longest to shortest until the current addition
#exceeds X/100 of the total contig length $contiglensum or custom threshold
my $Xstr = $customsum > 0 ? "G" . $X : $X;
$customsum = $customsum > 0 ? $customsum : $contiglensum;
print STDERR "Now calculating N", $Xstr, " value\n" if $debug;
my $curlen = 0;
my $curctg = 0; #For LX
OUTERLOOP: for (my $l = 0; $l < $arrlen; $l++) {
   for (my $m = 0; $m < $contiglens{$sortedcontiglens[$l]}; $m++) {
      $curctg++;
      $curlen = $curlen + $sortedcontiglens[$l];
      if ($curlen >= ($customsum * $X / 100)) {
         print "N", $Xstr, ": ", $sortedcontiglens[$l], "\n";
         print "L", $Xstr, ": ", $curctg, "\n";
         last OUTERLOOP;
      }
   }
}
#Extra features:
if ($numseqs > 0) {
   print "Shortest contig: ", $sortedcontiglens[$#sortedcontiglens], "\n";
   print "Longest contig: ", $sortedcontiglens[0], "\n";
   print "Number of contigs: ", $numseqs, "\n";
   print "Average contig: ", $contiglensum / $numseqs, "\n";
} else {
   print "N", $Xstr, ": 0\n";
   print "L", $Xstr, ": 0\n";
   print "Shortest contig: 0\n";
   print "Longest contig: 0\n";
   print "Number of contigs: 0\n";
   print "Average contig: 0\n";
}
