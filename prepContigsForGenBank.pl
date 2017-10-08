#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#########################################################################################
# prepContigsForGenBank.pl                                                              #
# Version 1.0 (2017/04/20)                                                              #
# Version 1.1 (2017/05/10)                                                              #
# Description:                                                                          #
# This script annotates FASTQ headers for GenBank submission, and produces an AGP       #
#                                                                                       #
# Usage:                                                                                #
#  prepContigsForGenBank.pl [-i input_contigs.fastq] <Configuration String>                    #
# Options:                                                                              #
#  --input_file,-i:   		Input contigs FASTQ file name (default: STDIN)          #
#  --agp_file,-a:   		Output AGP file name                                     #
#  --gap_size,-g:		Gap size to use for AGP gap lines                        #
#  --gap_type,-t:		Gap type to use (N for known size, U for unknown size)   #
#  --mtDNA,-m:			Annotate mtDNA with these fields                         #
#
#  Configuration String:	Format string composed as follows:                      #
#				[Scaffold name]:[contig name 1]->[contig name 2]->      #
#				[contig name 3]*->[contig name 4]=>[Scaffold name]      #
#########################################################################################

=pod

=head1 NAME

prepContigsForGenBank.pl - Generate annotated FASTQ and AGP for GenBank submission

=head1 SYNOPSIS

prepContigsForGenBank.pl [options] <Configuration String>

 Options:
  --help,-h,-?		Display this help documentation
  --input_file,-i	Input contigs FASTQ file name (default: STDIN)
  --agp_file,-a		Output AGP file
  --unscaffolded,-u	Include unscaffolded contigs in the AGP
  --gap_size,-g         Gap size to use for AGP gap lines (Default: 500)
  --gap_type,-t         Gap type to use, either Unknown size, or kNown size (Default: N)
  --mtDNA,-m:		Annotate mtDNA with these fields

 Mandatory:
  Configuration String	Format string composed as follows:
			[Scaffold name]:[contig name 1]->[contig name 2]->
			...->[contig name k]*->...->[contig name n]=>
			[Scaffold name]:[contig name 1]->[contig name 2]*
			Note: The asterisk denotes reverse complementing the
			contig.
			This is used in assigning to chromosomes, and
                        generating the AGP

=head1 DESCRIPTION
This script annotates FASTQ headers for GenBank submission and generates an AGP
file based on the configuration string.  The configuration string consists of contig
names with * indicating the reverse complement of a contig, and -> as the delimiter
between contigs in a scaffold.  Each scaffold is prefixed with the desired output
name and a colon (:), and scaffolds are delimited by =>.
Outputs the annotated FASTQ to STDOUT.

=cut

my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $unscaffolded = 0;
my $config_string = "";
my $agp_file = "";
my $gap_type = "N"; #For now, we'll use specified size gaps
my $gap_size = 500; #For now, we'll specify the gap size as 500 nt
#If $gap_type is "U", $gap_size must be 100
my $gap_reason = "scaffold"; #We aren't specifying centromeres, telomeres, or heterochromatin
my $linkage_yn = "yes";
my $linkage_evidence = "align_genus"; #For Drosophila, we used MUMmer
#Perhaps use "align_xgenus" for protein synteny-based linkage
#Use "map" for linkage map, e.g. MSG linkage LOD
my $mtDNA_fields = "";
GetOptions('input_file|i=s' => \$input_path, 'agp_file|a=s' => \$agp_file, 'unscaffolded|u' => \$unscaffolded, 'gap_size|g=i' => \$gap_size, 'gap_type|t=s' => \$gap_type, 'mtDNA|m=s' => \$mtDNA_fields, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

if ($input_path ne "STDIN") {
   unless(open(CONTIGS, "<", $input_path)) {
      print STDERR "Error opening input contigs FASTQ file.\n";
      exit 2;
   }
} else {
   open(CONTIGS, "<&", "STDIN"); #Duplicate the file handle for STDIN to CONTIGS so we can seamlessly handle piping
}

if (scalar @ARGV < 1) { #Not enough mandatory arguments
   print STDERR "Missing Configuration String.\n";
   exit 3;
} else {
   $config_string = $ARGV[0];
}

if ($agp_file eq "") {
   print STDERR "Missing AGP file path for output.\n";
   exit 5;
} else {
   unless(open(AGP, ">", $agp_file)) {
      print STDERR "Error opening output AGP file.\n";
      exit 4;
   }
   #Make the AGP header:
   print AGP "##agp-version", "\t", "2.0", "\n";
   #Later perhaps add the ORGANISM, TAX_ID, ASSEMBLY NAME, ASSEMBLY DATE, etc. headers
}

chomp $config_string;

#Read in all the contigs (temporary solution):
my %contigs = ();
my %contig_quals = ();
my ($header, $sequence, $qualheader, $quals) = ('', '', '', '');
my $fastq_line_modulus = 1;
while (my $line = <CONTIGS>) {
   chomp $line;
   if ($fastq_line_modulus == 1) {
      $header = substr $line, 1;
   } elsif ($fastq_line_modulus == 2) {
      $sequence = $line;
   } elsif ($fastq_line_modulus == 3) {
      $qualheader = substr $line, 1;
   } else { #Modulus is 0
      $quals = $line;
      $contigs{$header} = $sequence;
      $contig_quals{$header} = $quals;
      ($header, $sequence, $qualheader, $quals) = ('', '', '', '');
   }
   $fastq_line_modulus++;
   $fastq_line_modulus %= 4;
}
#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close(CONTIGS);
}
#Keep track of the unscaffolded contigs if asked to:
my %scaffolded_contigs = ();

#Keep track of the mtDNA contig if it exists:
my $mtDNA = "";
my ($mtDNA_prefix, $mtDNA_contig, $mtDNA_scaf, $mtDNA_fieldcodes);
if ($mtDNA_fields ne "") {
   ($mtDNA_prefix, $mtDNA_contig) = split /:/, $mtDNA_fields, 2;
   ($mtDNA_scaf, $mtDNA_fieldcodes) = split /\s+/, $mtDNA_prefix, 2;
}

#Decompose the configuration string:
my @scaffold_arr = split /=>/, $config_string; #Into scaffolds
for my $scaffold (@scaffold_arr) {
   my $scaffold_start = 1;
   my $scaffold_part_num = 1;
   my ($prefix, $contig_string) = split /:/, $scaffold; #Into prefix and contig string
   my @contig_arr = split /->/, $contig_string; #Into contigs
   for my $contig (@contig_arr) {
      my $contig_name = $contig;
      my $orientation = "?";
      if ($contig =~ /\*$/) {
         $orientation = "-";
         $contig_name = substr $contig_name, 0, -1;
      } else {
         $orientation = "+";
      }
      my $genbank_header = $contig_name; #Not sure if we need to add more here
      $genbank_header .= " ${mtDNA_fieldcodes}" if ($mtDNA_fields ne "" and $contig eq $mtDNA_contig);
      my $contig_start = 1;
      my $contig_end = length($contigs{$contig_name});
      my $agp_line = join("\t", $prefix, $scaffold_start, $scaffold_start+$contig_end-1, $scaffold_part_num, "W", $contig_name, $contig_start, $contig_end, $orientation) . "\n";
      $scaffold_start += $contig_end;
      $scaffold_part_num++;
      my $agp_gap_line = join("\t", $prefix, $scaffold_start, $scaffold_start+$gap_size-1, $scaffold_part_num, $gap_type, $gap_size, $gap_reason, $linkage_yn, $linkage_evidence) . "\n";
      $scaffold_start += $gap_size;
      $scaffold_part_num++;
      #Print the contig and gap lines to the AGP:
      print AGP $agp_line;
      print AGP $agp_gap_line unless $contig eq $contig_arr[-1];
      #Output the annotated FASTQ record to STDOUT:
      print STDOUT "@", $genbank_header, "\n", $contigs{$contig_name}, "\n", "+", $genbank_header, "\n", $contig_quals{$contig_name}, "\n";
      #Add the contig to the list of scaffolded contigs:
      $scaffolded_contigs{$contig_name} = 1;
   }
}

#Output the unscaffolded contigs if asked:
#Note that GenBank requires column 1 to be different from column 6, so
# simply append _scaf to the name for now.
if ($unscaffolded != 0) {
   for my $key (keys %contigs) {
      my $genbank_header = $key; #Not sure if we need to add more here
      #Output the AGP line for an unplaced contig:
      print AGP join("\t", $key . "_scaf", 1, length($contigs{$key}), 1, "W", $key, 1, length($contigs{$key}), "+"), "\n" unless exists($scaffolded_contigs{$key});
      #Output the annotated FASTQ record to STDOUT:
      print STDOUT "@", $genbank_header, "\n", $contigs{$key}, "\n", "+", $genbank_header, "\n", $contig_quals{$key}, "\n" unless exists($scaffolded_contigs{$key});
   }
}

close(AGP);

exit 0;
