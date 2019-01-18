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
# Version 1.0 (2018/11/21) Tested on Dyak NY73PB v2            #
################################################################

#First pass script to identify KEGG pathways for a given gene/transcript
# based on the output of Trinotate, since Trinotate occasionally
# only outputs a KEGG gene, but not the associated KO:K#####
#This script takes advantage of the KEGG REST API, so please don't spam it

my $SCRIPTNAME = "amendTrinotateKEGG.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

amendTrinotateKEGG.pl - Output KO terms for genes/transcripts using Trinotate report as input

=head1 SYNOPSIS

amendTrinotateKEGG.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_xls,-i         Path to input Trinotate report (default: STDIN)
  --genes,-g             Path to transcript output from this script
                         so that we output in terms of genes
                         (omit this if you want one line per transcript)
  --version,-v           Output version string

=head1 DESCRIPTION

This script outputs a list of KO terms for each gene/transcript, based
on the KEGG genes indicated by the input Trinotate report. It looks up
KO terms using the KEGG REST API, so please don't run this script too
often, or else the KEGG webmasters may get annoyed.
The output has a column for the gene or transcript ID, and the second
column is a comma-separated list of KO terms.

=cut

my $help = 0;
my $man = 0;
my $trinotate_path = "STDIN";
my $genes = "";
my $dispversion = 0;
GetOptions('input_xls|i=s' => \$trinotate_path, 'genes|g=s' => \$genes, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the trinotate report, or set it up to be read from STDIN:
my $trinotatefh;
if ($trinotate_path ne "STDIN") {
   unless(open($trinotatefh, "<", $trinotate_path)) {
      print STDERR "Error opening Trinotate report ${trinotate_path}.\n";
      exit 2;
   }
} else {
   open($trinotatefh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $trinotatefh so we can seamlessly handle piping
}

my $genesfh;
unless (-e $genes) {
   print STDERR "Input per-transcript KEGG file ${genes} does not exist, so generating per-transcript KEGG file on STDOUT.\n";
} else {
   unless(open($genesfh, "<", $genes)) {
      print STDERR "Error opening per-transcript KEGG file ${genes}.\n";
      close($trinotatefh);
      exit 2;
   }

}

my $ht = HTTP::Tiny->new;

my $header_line = <$trinotatefh>;
chomp $header_line;
my @header = split /\t/, $header_line;
my @kegg_indices = grep {$header[$_] eq "Kegg"} 0..$#header;
my $kegg_index = $kegg_indices[0];

#Iterate over each transcript in the Trinotate report, and compile a list of KO terms:
my %transcript_KOs = ();
my @input_transcripts = ();
my @input_genes = ();
my %gene_transcript_map = ();
while (my $line = <$trinotatefh>) {
   chomp $line;
   my @linearr = split /\t/, $line;
   my @KOs = ();
   #Skip the KEGG fetches if we're passed the genes file:
   unless (-e $genes) {
      #If the Trinotate report doesn't have any hits, leave the KO array empty:
      $linearr[$kegg_index] = "" if $linearr[$kegg_index] eq ".";
      #Otherwise, separate the elements of this column out:
      my @keggarr = split /`/, $linearr[$kegg_index];
      #Iterate over each element:
      for my $kegg_item (@keggarr) {
         #If the KO term is already annotated, store it in the KO array:
         if ($kegg_item =~ /^KO:/) {
            push @KOs, $kegg_item;
         } else { #Assuming that if not KO:, we are given a KEGG gene
            #Otherwise, if we are provided a KEGG gene, eliminate the prefix
            # and then look up the KO via the KEGG REST API:
            $kegg_item =~ s/KEGG://g;
            my $url = "http://rest.kegg.jp/get/" . $kegg_item;
            my $response = $ht->get($url);
            #Note: This only detects the first ORTHOLOGY line with a K prefix, so
            # in any cases with multiple KOs, only the first is retained.
            # I doubt there are such cases, but will still leave a note.
            if ($response->{success} and $response->{content} =~ /ORTHOLOGY\s+(K\d+)\s+/m) {
               push @KOs, "KO:" . $1;
            }
         }
      }
   }
   #These make sure we retain input gene and transcript order:
   push @input_genes, $linearr[0];
   push @input_transcripts, $linearr[1];
   #Keep a mapping of genes to transcripts in case we want to output per-gene:
   $gene_transcript_map{$linearr[0]} = [] unless exists($gene_transcript_map{$linearr[0]});
   push @{$gene_transcript_map{$linearr[0]}}, $linearr[1];
   #Store the KO IDs for this transcript if we fetched them:
   push @{$transcript_KOs{$linearr[1]}}, @KOs unless -e $genes;
}
close($trinotatefh);

#If we're passed a per-transcript KEGG file, process it:
if (-e $genes) {
   while (my $line = <$genesfh>) {
      chomp $line;
      my ($tx, $KO_list) = split /\t/, $line, 2;
      my @KOs = split /,/, $KO_list;
      push @{$transcript_KOs{$tx}}, @KOs;
   }
   close($genesfh);
}

if (-e $genes) {
   for my $gene (uniq(@input_genes)) {
      my @gene_KOs = ();
      for my $tx (@{$gene_transcript_map{$gene}}) {
         push @gene_KOs, @{$transcript_KOs{$tx}};
      }
      my @uniq_gene_KOs = uniq @gene_KOs;
      print join("\t", $gene, join(",", @uniq_gene_KOs)), "\n";
   }
} else {
   my @unique_transcripts = uniq @input_transcripts;
   for my $tx (@unique_transcripts) {
      my @gene_KOs = uniq @{$transcript_KOs{$tx}};
      print join("\t", $tx, join(",", @gene_KOs)), "\n";
   }
}

exit 0;
