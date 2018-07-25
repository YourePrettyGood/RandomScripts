#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#########################################################################################
# manualScaffold.pl                                                                     #
# Version 1.0 (2016/02/18)                                                              #
# Version 1.1 (2017/02/23)                                                              #
# Version 1.2 (2017/06/08) Added file option for config string                          #
# Version 1.3 (2018/07/25) Added proper interpretation of N and U in AGP                #
# Description:                                                                          #
# This script concatenates contigs together into a scaffold based on a configuration    #
# string, separating contigs by 500 Ns.  The configuration string consists of contig    #
# names with * indicating the reverse complement of a contig, and -> as the delimiter   #
# between contigs in a scaffold.  Each scaffold is prefixed with the desired output     #
# name and a colon (:), and scaffolds are delimited by =>.                              #
# As of version 1.3, if you supply an AGP, N and U gap records are interpreted according#
# to the AGP version 2.0 specification.                                                 #
#                                                                                       #
# Usage:                                                                                #
#  manualScaffold.pl [-i input_contigs.fasta] <Configuration String>                    #
# Options:                                                                              #
#  --input_file,-i:   		Input contigs FASTA file name (default: STDIN)          #
#  --config_string,-c:		File containing long configuration string               #
#  --agp_file,-a:   		Input AGP file name                                     #
#
#  Configuration String:	Format string composed as follows:                      #
#				[Scaffold name]:[contig name 1]->[contig name 2]->      #
#				[contig name 3]*->[contig name 4]=>[Scaffold name]      #
#########################################################################################

=pod

=head1 NAME

manualScaffold.pl - Generate scaffolds by manually joining contigs

=head1 SYNOPSIS

manualScaffold.pl [options] <Configuration String>

 Options:
  --help,-h,-?		Display this help documentation
  --input_file,-i	Input contigs FASTA file name (default: STDIN)
  --config_string,-c    File containing long configuration string (optional)
  --agp_file,-a		Input AGP file (instead of configuration string)
  --unscaffolded,-u	Output unscaffolded contigs individually at the end

 Mandatory:
  Configuration String	Format string composed as follows:
			[Scaffold name]:[contig name 1]->[contig name 2]->
			...->[contig name k]*->...->[contig name n]=>
			[Scaffold name]:[contig name 1]->[contig name 2]*
			Note: The asterisk denotes reverse complementing the
			contig.
			An AGP file may be used instead of the configuration
			string.

=head1 DESCRIPTION
This script concatenates contigs together into a scaffold based on a configuration
string, separating contigs by 500 Ns.  The configuration string consists of contig
names with * indicating the reverse complement of a contig, and -> as the delimiter
between contigs in a scaffold.  Each scaffold is prefixed with the desired output
name and a colon (:), and scaffolds are delimited by =>.
As of version 1.3, if you supply an AGP, this script will interpret N and U records
as per the AGP version 2.0 specification, rather than forcing N of 500 bp.
Outputs the scaffolds to STDOUT.

=cut

my $scriptname = "$0";
my $version = "1.3";

my $display_version = 0;
my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $unscaffolded = 0;
my $config_string = "";
my $agp_file = "";
my $config_string_file = "";
GetOptions('input_file|i=s' => \$input_path, 'config_string|c=s' => \$config_string_file, 'agp_file|a=s' => \$agp_file, 'unscaffolded|u' => \$unscaffolded, 'help|h|?' => \$help, man => \$man, 'version|v' => \$display_version) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;
die join(" ", $scriptname, "version", $version) if $display_version;

if ($input_path ne "STDIN") {
   unless(open(CONTIGS, "<", $input_path)) {
      print STDERR "Error opening input contigs FASTA file.\n";
      exit 2;
   }
} else {
   open(CONTIGS, "<&", "STDIN"); #Duplicate the file handle for STDIN to CONTIGS so we can seamlessly handle piping
}

my %scaffold_gaps = (); #Store gaps between contigs from AGP
my $default_gap_size = 500;

if ($agp_file eq "") {
   if (scalar @ARGV < 1) { #Not enough mandatory arguments
      if ($config_string_file ne "") {
         unless(open(CONFIGSTR, "<", $config_string_file)) {
            print STDERR "Error opening long configuration string file ${config_string_file}.\n";
            exit 4;
         }
         $config_string = <CONFIGSTR>;
      } else {
         print STDERR "Missing Configuration String.\n";
         exit 3;
      }
   } else {
      $config_string = $ARGV[0];
   }
} else {
   $default_gap_size = 100; #AGP 2.0 says U records get 100 bp gaps
   unless(open(AGP, "<", $agp_file)) {
      print STDERR "Error opening input AGP file.\n";
      exit 4;
   }
   #Convert the AGP file to a configuration string:
   my @chroms = ();
   my %confighash = ();
   my $prev_contig = "";
   while (my $line = <AGP>) {
      chomp $line;
      next if $line =~ /^#/;
      my @line_elements = split /\t/, $line;
      next if $line_elements[4] !~ /[DNUW]/;
      if (scalar@line_elements < 9) {
         print STDERR "Bad D, N, U, or W record ", $line, "\n";
         next;
      }
      if ($line_elements[4] =~ /[DW]/) {
         push @chroms, $line_elements[0] if scalar@chroms == 0 or $chroms[$#chroms] ne $line_elements[0];
         my $config_part = $line_elements[5] . ($line_elements[8] eq "-" ? "*" : "");
         if (exists($confighash{$line_elements[0]})) {
            $confighash{$line_elements[0]} .= "->";
         } else {
            $confighash{$line_elements[0]} = $line_elements[0] . ":";
         }
         $confighash{$line_elements[0]} .= $config_part;
         $prev_contig = $config_part;
      } else { #Gap record
         if ($line_elements[4] eq "N") { #Specified gap size
            $scaffold_gaps{$prev_contig} = $line_elements[5]+0; #Make sure it's numeric
         } else { #U record, so unknown gap size, default is 100 bp
            $scaffold_gaps{$prev_contig} = $default_gap_size;
         }
      }
   }
   my $output = "";
   for my $chrom (@chroms) {
      $output .= "=>" unless $output eq "";
      $output .= $confighash{$chrom};
   }
   $config_string = $output;
   close(AGP);
}

chomp $config_string;

#Read in all the contigs (temporary solution):
my %contigs = ();
my ($header, $sequence) = ('', '');
while (my $line = <CONTIGS>) {
   chomp $line;
   #Put the contig in the hash if we've reached the next contig:
   if ($header ne '' and $line =~ />/) {
      $contigs{$header} = $sequence;
      ($header, $sequence) = ('', '');
   }
   #Fill up the header and sequence variables if we haven't reached the next contig:
   if ($line =~ />/) {
      $header = substr $line, 1;
      if ($line =~ /\s+/) { #If header has any spaces, only use the first word
         my @header_words = split /\s+/, $header;
         $header = $header_words[0];
      }
   } else {
      $sequence .= $line;
   }
}
#Add the final contig into the hash, if it exists:
if ($header ne '' and $sequence ne '') {
   $contigs{$header} = $sequence;
   ($header, $sequence) = ('', '');
}
#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close(CONTIGS);
}
#Keep track of the unscaffolded contigs if asked to:
my %scaffolded_contigs = ();

#Decompose the configuration string:
my @scaffold_arr = split /=>/, $config_string; #Into scaffolds
for my $scaffold (@scaffold_arr) {
   my ($prefix, $contig_string) = split /:/, $scaffold; #Into prefix and contig string
   my @contig_arr = split /->/, $contig_string; #Into contigs
   my $scaffold_sequence = ''; #Build the scaffold sequence up in this variable
   for my $contig (@contig_arr) {
      #If no AGP provided, add 500 N gaps between contigs:
      $scaffold_sequence .= "N"x$default_gap_size if $scaffold_sequence ne "" and $agp_file eq "";
      if ($contig =~ /\*$/) {
         #Reverse complement the sequence, and save the revcomp:
         my $hashkey = substr $contig, 0, -1;
         $scaffolded_contigs{$hashkey} = 1;
         my $complement = $contigs{$hashkey};
         $complement =~ tr/ACGTNacgtn/TGCANtgcan/;
         $scaffold_sequence .= reverse $complement;
      } else {
         #Save the sequence:
         $scaffolded_contigs{$contig} = 1;
         $scaffold_sequence .= $contigs{$contig};
      }
      #Add the N spacer after the contig if an AGP was provided:
      $scaffold_sequence .= "N"x$scaffold_gaps{$contig} if exists($scaffold_gaps{$contig});
   }
   print STDOUT ">${prefix}\n${scaffold_sequence}\n";
}

#Output the unscaffolded contigs if asked:
if ($unscaffolded != 0) {
   for my $key (keys %contigs) {
      print STDOUT ">${key}\n", $contigs{$key}, "\n" unless exists($scaffolded_contigs{$key});
   }
}

exit 0;
