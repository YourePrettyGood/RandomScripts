#!/usr/bin/env perl

use warnings;
use strict;
use HTTP::Tiny;
use List::MoreUtils qw(uniq);
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.0 (2018/11/27) Tested on Dyak NY73PB v2            #
################################################################

#First pass script to retrieve KEGG pathway descriptions based on
# the output of amendTrinotateKEGG.pl
#This script takes advantage of the KEGG REST API, so please don't spam it

my $SCRIPTNAME = "getKEGGdescriptions.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

getKEGGdescriptions.pl - Retrieve KEGG pathway descriptions for KOs from amendTrinotateKEGG.pl

=head1 SYNOPSIS

amendTrinotateKEGG.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input,-i             Path to input gene KO mapping file (default: STDIN)
  --version,-v           Output version string
  --debug,-d             Output debugging info to STDERR

=head1 DESCRIPTION

This script outputs a TSV of KO pathways and their descriptions, as
retrieved using the KEGG REST API. The expected input is the output
of amendTrinotateKEGG.pl, which consists of two tab-separated columns,
the first being the gene or transcript ID, and the second being a
comma-separated list of KEGG pathway IDs (KO:#####).
The output is simply the KEGG pathway ID, and its description.

=cut

my $help = 0;
my $man = 0;
my $kegg_path = "STDIN";
my $dispversion = 0;
my $debug = 0;
GetOptions('input|i=s' => \$kegg_path, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the gene-KEGG mapping file, or set it up to be read from STDIN:
my $keggfh;
if ($kegg_path ne "STDIN") {
   unless(open($keggfh, "<", $kegg_path)) {
      print STDERR "Error opening gene-KEGG map ${kegg_path}.\n";
      exit 2;
   }
} else {
   open($keggfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $keggfh so we can seamlessly handle piping
}

my $ht = HTTP::Tiny->new;

#Read in all the KEGG pathway IDs into a hash so we eliminate redundancy:
my %KOs = ();
print STDERR "Reading gene-KEGG mapping file.\n";
while (my $line = <$keggfh>) {
   chomp $line;
   my ($gene, $KOlist) = split /\t/, $line, 2;
   next if $KOlist !~ /KO:/i;
   my @KO_arr = split /,/, $KOlist;
   for my $KO (@KO_arr) {
      $KOs{$KO} = 1;
   }
}
close($keggfh);
print STDERR "Done reading gene-KEGG mapping file.\n";

print STDERR "Retrieving descriptions for ", scalar(keys %KOs), " KEGG pathways.\n";
for my $KO (keys %KOs) {
   my ($prefix, $kid) = split /:/, $KO, 2;
   my $url = "http://rest.kegg.jp/get/" . lc($prefix) . ":" . $kid;
   my $response = $ht->get($url);
   #Note: This only detects the first ORTHOLOGY line with a K prefix, so
   # in any cases with multiple KOs, only the first is retained.
   # I doubt there are such cases, but will still leave a note.
   print STDERR join("\t", $KO, $1), "\n" if $debug and $response->{success} and $response->{content} =~ /DEFINITION\s+(.+)/m;
   print STDERR $url, "\n" if $debug > 2;
   print STDERR $response->{content}, "\n" if $debug > 2;
   print STDERR $response->{content}, "\n" if $debug > 1 and $response->{success};
   if ($response->{success} and $response->{content} =~ /^DEFINITION\s+(.+)$/m) {
      print STDOUT join("\t", $KO, $1), "\n";
   }
}
print STDERR "Done retrieving descriptions for KEGG pathways.\n";

exit 0;
