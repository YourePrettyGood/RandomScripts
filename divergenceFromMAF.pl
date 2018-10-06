#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to calculate divergence from pairwise MAF
# of 1:1 alignments produced by LAST

my $SCRIPTNAME = "divergenceFromMAF.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

divergenceFromMAF.pl - Calculate divergence from pairwise 1:1 MAF

=head1 SYNOPSIS

divergenceFromMAF.pl [options] <Species 1 prefix> <Species 2 prefix>

 Options:
  --help,-h,-?          Print this help documentation
  --input_MAF,-i        Path to input MAF 1:1 alignment file (default: STDIN)
  --sort,-s             Sort output by scaffold and position
  --debug,-d            Output debugging information to STDERR
  --version,-v          Output version string

=head1 DESCRIPTION

This script calculates divergence between a pair of species for
sites covered by 1:1 alignments in MAF format. Coordinates of
sites are in the space of the first species. The scaffold IDs in
the MAF file should follow the format [species prefix].[scaffold name],
where the species prefix and scaffold name both are alphanumeric
(plus underscores). As long as the scaffold IDs follow this convention,
and you provide the species prefixes precisely as found in the MAF,
the script will work as intended.

For example:

=begin text

a score=10 mismap=1e-05
s DyakTai18E2.2L 300 10 + 29953808 AATAACGGCT
s DyakNY73PB.2L  300 10 + 24234981 AACAACGGGT
p                                  !!!!!!!!!!
p                                  !!!!!!!!!!

=end text

should output:

=begin text

2L	301	0	0
2L	302	0	0
2L	303	1	0
2L	304	0	0
2L	305	0	0
2L	306	0	0
2L	307	0	0
2L	308	0	0
2L	309	1	0
2L	310	0	0

=end text

Sites not covered by alignments will not be output, so you
will have to post-process the output to add lines with 4th
column equal to 1 for unaligned sites.

=cut

my $help = 0;
my $man = 0;
my $maf_path = "STDIN";
my $sort_output = 0;
my $debug = 0;
my $dispversion = 0;
GetOptions('input_maf|i=s' => \$maf_path, 'sort|s' => \$sort_output, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the MAF file, or set it up to be read from STDIN:
print STDERR "Opening MAF file\n" if $debug;
my $maf_fh;
if ($maf_path ne "STDIN") {
   unless(open($maf_fh, "<", $maf_path)) {
      print STDERR "Error opening MAF file: ${maf_path}.\n";
      exit 2;
   }
} else {
   open($maf_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $maf_fh so we can seamlessly handle piping
}

#Extract the positional arguments, which are prefixes for the species:
my @prefixes = @ARGV;
#Make a hash so we can quickly check if the species from an "s" line matches:
print STDERR "Constructing prefix hash\n" if $debug;
my %prefix_hash = map {$_ => 1} @prefixes;

if (scalar(@prefixes) != 2) {
   close($maf_fh);
   print STDERR "It appears you have provided more than 2 prefixes. We only support pairwise alignments for the moment.\n";
   exit 3;
}

#Hash of hashes to store the elements of a MAF alignment:
my %aligned_seqs = ();
#Hash of offsets to stay in the appropriate coordinate space:
my %offsets = ();
#Array of species in the alignment:
my @species_aligned = ();
#Number of p lines we've encountered, which should equal the number of s lines:
my $num_p_lines = 0;

#Keep track of the total number of alignments for debugging:
my $num_a_lines = 0;

print STDERR "Parsing MAF ${maf_path}\n" if $debug;
my %output_blocks = ();
my @block_keys = ();
while (my $line = <$maf_fh>) {
   chomp $line;
   next if $line =~ /^#/ or $line eq ""; #Skip header and empty lines
   my @maf_elems = split /\s+/, $line; #MAF uses padded spacing
   @species_aligned = () if $maf_elems[0] eq "a"; #Reset species on new alignment
   $num_p_lines = 0 if $maf_elems[0] eq "a"; #Reset p count on new alignment
   $num_a_lines++ if $maf_elems[0] eq "a";
   print STDERR "Parsed ${num_a_lines} alignments\n" if $debug and ${num_a_lines} % 1000 == 1 and $maf_elems[0] eq "a";
   if ($maf_elems[0] eq "s") { #For aligned Sequence records, store the details
      my ($species, $scaffold) = split /\./, $maf_elems[1], 2;
      push @species_aligned, $species;
      my @seq_arr = split //, uc($maf_elems[6]);
      $aligned_seqs{$species} = {'scaffold' => $scaffold,
         'start' => $maf_elems[4] eq "-" ? $maf_elems[5]-$maf_elems[2] : $maf_elems[2]+1,
         'strand' => $maf_elems[4],
         'seq' => \@seq_arr};
   }
   if ($maf_elems[0] eq "p") { #Identify SNPs once we're in the Probability lines
      $num_p_lines++;
      next unless $num_p_lines == scalar(@species_aligned); #Skip unless we're on the last p line
      #Make a hash of the aligned species so we can take set differences:
      my %aligned_spp = map {$_ => 1} @species_aligned;

      #Sequence length for each "s" record should be identical, so just take
      # the length from the first species in the alignment:
      my $aln_length = scalar(@{$aligned_seqs{$species_aligned[0]}{'seq'}});

      #Initialize the offsets for each coordinate space:
      for my $species (@species_aligned) {
         $offsets{$species} = 0;
      }

      my $block_key = join(":", $aligned_seqs{$prefixes[0]}{'scaffold'}, $aligned_seqs{$prefixes[0]}{'start'});
      $output_blocks{${block_key}} = "";
      #Iterate over each column of the alignment to ID SNPs:
      for (my $i = 0; $i < $aln_length; $i++) {
         #ID divergent or shared sites:
         my $divergent = 0;
         my $omit_position = 0;
         #Skip deletions in species 1:
         next if $aligned_seqs{$prefixes[0]}{'seq'}[$i] eq '-';
         #Omit position if insertion in species 1:
         $omit_position = 1 if $aligned_seqs{$prefixes[1]}{'seq'}[$i] eq '-';
         $divergent = 1 if $aligned_seqs{$prefixes[0]}{'seq'}[$i] ne $aligned_seqs{$prefixes[1]}{'seq'}[$i];
         #Store the state of the position:
         my $species1_position = $aligned_seqs{$prefixes[0]}{'strand'} eq "-" ? $aligned_seqs{$prefixes[0]}{'start'}-$offsets{$prefixes[0]} : $aligned_seqs{$prefixes[0]}{'start'}+$offsets{$prefixes[0]};
         if ($aligned_seqs{$prefixes[0]}{'strand'} eq "+") {
            $output_blocks{${block_key}} .= join("\t", $aligned_seqs{$prefixes[0]}{'scaffold'}, $species1_position, $divergent, $omit_position) . "\n";
         } else {
            $output_blocks{${block_key}} = join("\t", $aligned_seqs{$prefixes[0]}{'scaffold'}, $species1_position, $divergent, $omit_position) . "\n" . $output_blocks{${block_key}};
         }
         for my $species (@species_aligned) {
            $offsets{$species}++ unless $aligned_seqs{$species}{'seq'}[$i] eq "-";
         }
      }
      push @block_keys, $block_key if $sort_output;
      print STDOUT $output_blocks{${block_key}} unless $sort_output;
   }
}
close($maf_fh);
print STDERR "Done parsing all ${num_a_lines} alignments from MAF file\n" if $debug;

if ($sort_output) {
   print STDERR "Sorting output blocks\n" if $debug;
   my @sorted_block_keys = map { join(":", @{$_}) } sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } map { [ split(":", $_) ] } @block_keys;
   print STDERR "Printing output blocks\n" if $debug;
   for my $block_key (@sorted_block_keys) {
      print STDOUT $output_blocks{${block_key}};
   }
}
print STDERR "Done printing output blocks\n" if $debug;

exit 0;
