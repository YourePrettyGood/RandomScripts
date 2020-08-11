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
# Version 1.4 (2018/11/15) Handle of IUPAC bases in revcomp, and keep input contig order#
# Version 1.5 (2019/10/20) Allow customization of fixed gap size for config string mode #
#                          Enables GenBank-compatible scaffolding (gap_size=100)        #
# Description:                                                                          #
# This script concatenates contigs together into a scaffold based on a configuration    #
# string, separating contigs by 500 Ns.  The configuration string consists of contig    #
# names with * indicating the reverse complement of a contig, and -> as the delimiter   #
# between contigs in a scaffold.  Each scaffold is prefixed with the desired output     #
# name and a colon (:), and scaffolds are delimited by =>.                              #
# As of version 1.3, if you supply an AGP, N and U gap records are interpreted according#
# to the AGP version 2.0 specification.                                                 #
# Version 1.5 enables modifying the fixed gap size in config string mode, so that       #
# -g 100 would fit with GenBank expectations.                                           #
#                                                                                       #
# Usage:                                                                                #
#  manualScaffold.pl [-i input_contigs.fasta] <Configuration String>                    #
# Options:                                                                              #
#  --help,-h,-?          Display this help documentation                                 #
#  --input_file,-i:      Input contigs FASTA file name (default: STDIN)                  #
#  --config_string,-c:   File containing long configuration string                       #
#  --agp_file,-a:        Input AGP file name                                             #
#  --unscaffolded,-u	 Output unscaffolded contigs individually at the end             #
#  --gap_size,-g         Fixed gap size to insert (for config string only)              #
#                        (Default: 500)                                                  #
#  --version,-v          Output version string                                          #
#                                                                                       #
#  Configuration String:	Format string composed as follows:                      #
#				[Scaffold name]:[contig name 1]->[contig name 2]->      #
#				[contig name 3]*->[contig name 4]=>[Scaffold name]      #
#########################################################################################

my $SCRIPTNAME = "manualScaffold.pl";
my $VERSION = "1.5";

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
  --gap_size,-g         Fixed gap size to insert (for config string only)
                        (Default: 500)
  --version,-v          Output version string

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
Version 1.5 adds the ability to specify a fixed gap size other than 500 when using
a config string. However, AGP mode still uses 100 Ns for U records.
Outputs the scaffolds to STDOUT.

=cut

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/AaCcGgTtRrYySsWwKkMmBbDdHhVvNn/TtGgCcAaYyRrSsWwMmKkVvHhDdBbNn/; #Complement incl. IUPAC degenerate bases
   return $reverse_sequence;
}

my $display_version = 0;
my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $unscaffolded = 0;
my $config_string = "";
my $agp_file = "";
my $config_string_file = "";
my $default_gap_size = 500;
my $dispversion = 0;
GetOptions('input_file|i=s' => \$input_path, 'config_string|c=s' => \$config_string_file, 'agp_file|a=s' => \$agp_file, 'unscaffolded|u' => \$unscaffolded, 'gap_size|g=i' => \$default_gap_size, 'help|h|?+' => \$help, man => \$man, 'version|v' => \$dispversion) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

print STDERR "Invalid gap size ${default_gap_size}, must be non-negative.\n" unless $default_gap_size >= 0;
exit 5 unless $default_gap_size >= 0;

my $contigsfh;
if ($input_path ne "STDIN") {
   unless(open($contigsfh, "<", $input_path)) {
      print STDERR "Error opening input contigs FASTA file.\n";
      exit 2;
   }
} else {
   open($contigsfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $contigsfh so we can seamlessly handle piping
}

my %scaffold_gaps = (); #Store gaps between contigs from AGP

if ($agp_file eq "") {
   if (scalar @ARGV < 1) { #Not enough mandatory arguments
      if ($config_string_file ne "") {
         my $configstrfh;
         unless(open($configstrfh, "<", $config_string_file)) {
            print STDERR "Error opening long configuration string file ${config_string_file}.\n";
            exit 4;
         }
         $config_string = <$configstrfh>;
         close($configstrfh);
      } else {
         print STDERR "Missing Configuration String.\n";
         exit 3;
      }
   } else {
      $config_string = $ARGV[0];
   }
} else {
   $default_gap_size = 100; #AGP 2.0 says U records get 100 bp gaps
   my $agpfh;
   unless(open($agpfh, "<", $agp_file)) {
      print STDERR "Error opening input AGP file.\n";
      exit 4;
   }
   #Convert the AGP file to a configuration string:
   my @chroms = ();
   my %confighash = ();
   my $prev_contig = "";
   while (my $line = <$agpfh>) {
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
   close($agpfh);
}

chomp $config_string;

#Read in all the contigs (temporary solution):
my @contigorder = ();
my %contigs = ();
my ($header, $sequence) = ('', '');
while (my $line = <$contigsfh>) {
   chomp $line;
   #Put the contig in the hash if we've reached the next contig:
   if ($header ne '' and $line =~ />/) {
      $contigs{$header} = $sequence;
      push @contigorder, $header;
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
   push @contigorder, $header;
   ($header, $sequence) = ('', '');
}
#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close($contigsfh);
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
         $scaffold_sequence .= revcomp($contigs{$hashkey});
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

my @unscaffolded_order = grep { !exists($scaffolded_contigs{$_}) } @contigorder;

#Output the unscaffolded contigs if asked:
if ($unscaffolded != 0) {
   for my $key (@unscaffolded_order) {
      print STDOUT ">${key}\n", $contigs{$key}, "\n";
   }
}

exit 0;
