#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

my $SCRIPTNAME = "configStringToAGP.pl";
my $VERSION = "1.2";

#Version 1.2 adds a flag for GenBank compatibility:
# Appends _scaf to the scaffold ID (column 1) of unscaffolded contigs

=pod

=head1 NAME

configStringToAGP.pl - Convert a config string for manualScaffold.pl to AGP

=head1 SYNOPSIS

configStringToAGP.pl [options] <Configuration String>

 Options:
  --help,-h,-?		Display this help documentation
  --input_fai,-i	Input contigs FASTA index file name (default: STDIN)
  --config_string,-c    File containing long configuration string (optional)
  --agp_file,-a		Output AGP file
  --unscaffolded,-u	Output unscaffolded contigs individually at the end
  --genbank_scaf,-s     Append _scaf to unscaffolded contigs (for GenBank)
  --gap_size,-g         Fixed gap size used from the config string
                        (default: 500)

 Mandatory:
  Configuration String	Format string composed as follows:
			[Scaffold name]:[contig name 1]->[contig name 2]->
			...->[contig name k]*->...->[contig name n]=>
			[Scaffold name]:[contig name 1]->[contig name 2]*
			Note: The asterisk denotes reverse complementing the
			contig.

=head1 DESCRIPTION
This script converts a configuration string as used by manualScaffold.pl
to AGP for easier visualization and identification of scaffolding
boundaries.

=cut

my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $unscaffolded = 0;
my $config_string = "";
my $agp_file = "STDOUT";
my $config_string_file = "";
my $gap_size = 500;
my $genbank_scaf = 0;
my $dispversion = 0;
GetOptions('input_file|i=s' => \$input_path, 'config_string|c=s' => \$config_string_file, 'agp_file|a=s' => \$agp_file, 'unscaffolded|u' => \$unscaffolded, 'gap_size|g=i' => \$gap_size, 'genbank_scaf|s' => \$genbank_scaf, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

print STDERR "Invalid gap size ${gap_size}, must be non-negative.\n" unless $gap_size >= 0;
exit 5 unless $gap_size >= 0;

my $contigs;
if ($input_path ne "STDIN") {
   unless(open($contigs, "<", $input_path)) {
      print STDERR "Error opening input contigs FASTA file.\n";
      exit 2;
   }
} else {
   open($contigs, "<&", "STDIN"); #Duplicate the file handle for STDIN to CONTIGS so we can seamlessly handle piping
}
my $agp;
if ($agp_file ne "STDOUT") {
   unless(open($agp, ">", $agp_file)) {
      print STDERR "Error opening output AGP file.\n";
      exit 3;
   }
} else {
   open($agp, ">&", "STDOUT"); #Duplicate the file handle for STDOUT to AGP so we can seamlessly handle piping
}

if (scalar @ARGV < 1) { #Not enough mandatory arguments
   if ($config_string_file ne "") {
      my $configstr;
      unless(open($configstr, "<", $config_string_file)) {
         print STDERR "Error opening long configuration string file ${config_string_file}.\n";
         exit 4;
      }
      $config_string = <$configstr>;
   } else {
      print STDERR "Missing Configuration String.\n";
      exit 3;
   }
} else {
   $config_string = $ARGV[0];
}

chomp $config_string;

#Read in all the contigs lengths from the FAI file:
my %ctg_lengths = ();
while (my $line = <$contigs>) {
   chomp $line;
   my @linearr = split /\t/, $line;
   $ctg_lengths{$linearr[0]} = $linearr[1];
}
close($contigs);

print $agp "##agp-version\t2.0\n";
#Keep track of the unscaffolded contigs if asked to:
my %scaffolded_contigs = ();

#Decompose the configuration string:
my @scaffold_arr = split /=>/, $config_string; #Into scaffolds
for my $scaffold (@scaffold_arr) {
   my ($prefix, $contig_string) = split /:/, $scaffold; #Into prefix and contig string
   my @contig_arr = split /->/, $contig_string; #Into contigs
   my $partnum = 1;
   my $scafstart = 1;
   for my $contig (@contig_arr) {
      my $orient = "+";
      my $hashkey = $contig;
      if ($contig =~ /\*$/) {
         $hashkey = substr $contig, 0, -1;
         $orient = "-";
      }
      $scaffolded_contigs{$hashkey} = 1;
      print STDERR "Missing ${hashkey} from contig lengths\n" unless exists($ctg_lengths{$hashkey});
      print $agp join("\t", $prefix, $scafstart, $scafstart+$ctg_lengths{$hashkey}-1, $partnum, "W", $hashkey, "1", $ctg_lengths{$hashkey}, $orient), "\n";
      $scafstart += $ctg_lengths{$hashkey};
      $partnum++;
      unless ($contig eq $contig_arr[$#contig_arr]) { #No gap after last contig
         print $agp join("\t", $prefix, $scafstart, $scafstart+$gap_size-1, $partnum, "U", $gap_size, "scaffold", "yes", "align_genus"), "\n";
         $scafstart += $gap_size;
         $partnum++;
      }
   }
}

#Output the unscaffolded contigs if asked:
if ($unscaffolded != 0) {
   for my $key (keys %ctg_lengths) {
      my $scaf_ID = $key;
      my @scaf_ID_parts = split /\s+/, $scaf_ID;
      $scaf_ID_parts[0] .= "_scaf" if $genbank_scaf;
      print $agp join("\t", join(" ", @scaf_ID_parts), "1", $ctg_lengths{$key}, "1", "W", $key, "1", $ctg_lengths{$key}, "+"), "\n" unless exists($scaffolded_contigs{$key});
   }
}
close($agp);
exit 0;
