#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#This script takes a GFF3 in one coordinate space, plus an AGP, and outputs a
# GFF3 in the other coordinate space indicated by the AGP.
#For instance, if an AGP from scaffolds to chromosomes is provided, a
# scaffold-space GFF3 could be converted to a chromosome-space GFF3, or vice
# versa (with the --invert flag).
#If an AGP from contigs to scaffolds is provided, a contig-space GFF3 could be
# converted to a scaffold-space GFF3, or vice versa (with the --invert flag).
#Version 1.1 adds compatibility with converting coordinates for BED files

my $SCRIPTNAME = "ScaffoldGFF3.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

ScaffoldGFF3.pl - Convert GFF3 between coordinate spaces defined by an AGP

=head1 SYNOPSIS

ScaffoldGFF3.pl -a [AGP file] -g [Source GFF3 file] [-m <BED file>] [--invert] [--broken <Broken GFF3 file>] > [Destination GFF3 file]

 Options:
  --help,-h,-?	        Display this help documentation
  --version,-v          Output version string
  --invert,-i           Convert coordinates from the more-joined space to less-joined
  --broken,-b           GFF3 file of features broken during conversion to output
 Mandatory:
  --agp,-a		AGP file for e.g. mapping scaffolds to chromosomes
  --gff,-g		GFF3 file of features in e.g. scaffold-space
  --bed,-m              BED file of features in e.g. scaffold-space
                        (i.e. convert BED rather than GFF)

=head1 DESCRIPTION

ScaffoldGFF3.pl converts a GFF3 from one coordinate space into another
using an AGP file. This allows for easy conversion of annotations between
coordinate spaces (as long as only scaffolding has occurred).

=cut

#Initialize the input parameters
my $help = 0;
my $man = 0;
my $debug = 0;
my $input_agp = "";
my $input_gff = "";
my $input_bed = "";
my $invert = 0;
my $broken_gff3 = "";
my $dispversion = 0;

#Fetch the command line parameters
GetOptions('version|v' => \$dispversion, 'help|h|?+' => \$help, 'agp|a=s' => \$input_agp, 'gff|g=s' => \$input_gff, 'invert|i' => \$invert, 'broken|b=s' => \$broken_gff3, 'bed|m=s' => \$input_bed, 'debug|d+' => \$debug, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

pod2usage(-exitval => 2, -output => \*STDERR) if $input_agp eq "" or ($input_gff eq "" and $input_bed eq "");

#Make sure the required parameters are filled out correctly
unless (-e $input_agp) {
   print STDERR "The AGP input file ${input_agp} does not exist.\n";
   exit 3;
}
unless (($input_gff ne "" and -e $input_gff) or ($input_bed ne "" and -e $input_bed)) {
   print STDERR "The GFF3 (${input_gff}) or BED (${input_bed}) input file does not exist.\n";
   exit 4;
}

#Functions for converting between coordinates:
sub contigToScaffold($$$$$$) {
   my $feature_coordinates = shift @_;
   my @ctg_list = @{shift @_};
   my @ctg_intervals = @{shift @_};
   my @scaf_intervals = @{shift @_};
   my @ctg_orientations = @{shift @_};
   my $debug = shift @_;
   
   #Split the coordinate string into components:
   my ($feature_ctg, $feature_interval_strand) = split /:/, $feature_coordinates, 2;
   my ($feature_interval, $strand) = split /\(/, $feature_interval_strand, 2;
   $strand = substr($strand, 0, -1);
   my ($feature_start, $feature_end) = split /-/, $feature_interval, 2;
   
   #Find the appropriate index for matching intervals:
   my $ctg_index = $#ctg_intervals;
   my ($ctg_start, $ctg_end);
   while ($ctg_index >= 0) {
      ($ctg_start, $ctg_end) = split /-/, $ctg_intervals[$ctg_index], 2;
      last if $feature_ctg eq $ctg_list[$ctg_index] and $feature_start >= $ctg_start and $feature_end <= $ctg_end;
      #Boundary cases:
      print STDERR join("\t", $feature_start, $feature_end, $ctg_start, $ctg_end), "\n" if $debug > 1 and $feature_ctg eq $ctg_list[$ctg_index] and $feature_start < $ctg_start and $feature_end >= $ctg_start and $feature_end <= $ctg_end;
      return undef if $feature_ctg eq $ctg_list[$ctg_index] and $feature_start < $ctg_start and $feature_end >= $ctg_start and $feature_end <= $ctg_end;
      print STDERR join("\t", $feature_start, $feature_end, $ctg_start, $ctg_end), "\n" if $debug > 1 and $feature_ctg eq $ctg_list[$ctg_index] and $feature_start >= $ctg_start and $feature_start <= $ctg_end and $feature_end > $ctg_end;
      return undef if $feature_ctg eq $ctg_list[$ctg_index] and $feature_start >= $ctg_start and $feature_start <= $ctg_end and $feature_end > $ctg_end;
      $ctg_index--;
   }
   my ($scaf_part_start, $scaf_part_end) = split /-/, $scaf_intervals[$ctg_index], 2;
   
   #Convert the coordinates:
   my ($feature_start_scaf, $feature_end_scaf);
   if ($ctg_orientations[$ctg_index] eq "-") {
      $feature_start_scaf = $ctg_end - $feature_end + $scaf_part_start;
      $feature_end_scaf = $ctg_end - $feature_start + $scaf_part_start;
      $strand = $strand eq "-" ? "+" : "-";
   } else {
      #Simply adjust for shift of contig relative to scaffold:
      $feature_start_scaf = $feature_start - $ctg_start + $scaf_part_start;
      $feature_end_scaf = $feature_end - $ctg_start + $scaf_part_start;
   }
   print STDERR join("\t", $feature_start, $feature_end, $feature_start_scaf, $feature_end_scaf, $ctg_start, $ctg_end, $scaf_part_start, $scaf_part_end, $ctg_orientations[$ctg_index], $ctg_index), "\n" if $debug > 2;
   return "${feature_start_scaf}-${feature_end_scaf}(${strand})";
}

sub scaffoldToContig($$$$$$) {
   my $scaf_coordinates = shift @_;
   my @ctg_list = @{shift @_};
   my @ctg_intervals = @{shift @_};
   my @scaf_intervals = @{shift @_};
   my @ctg_orientations = @{shift @_};
   my $debug = shift @_;
   
   #Split the coordinate string into components:
   my ($feature_interval, $strand) = split /\(/, $scaf_coordinates, 2;
   $strand = substr($strand, 0, -1);
   my ($feature_start, $feature_end) = split /-/, $feature_interval, 2;
   
   #Find the corresponding contig:
   my $ctg_index = $#ctg_intervals;
   my ($scaf_part_start, $scaf_part_end);
   while ($ctg_index >= 0) {
      ($scaf_part_start, $scaf_part_end) = split /-/, $scaf_intervals[$ctg_index], 2;
      last if $feature_start >= $scaf_part_start and $feature_end <= $scaf_part_end;
      #Boundary cases:
      print STDERR join("\t", $feature_start, $feature_end, $scaf_part_start, $scaf_part_end), "\n" if $debug > 1 and $feature_start < $scaf_part_start and $feature_end >= $scaf_part_start and $feature_end <= $scaf_part_end;
      return undef if $feature_start < $scaf_part_start and $feature_end >= $scaf_part_start and $feature_end <= $scaf_part_end;
      print STDERR join("\t", $feature_start, $feature_end, $scaf_part_start, $scaf_part_end), "\n" if $debug > 1 and $feature_start >= $scaf_part_start and $feature_end <= $scaf_part_start and $feature_end > $scaf_part_end;
      return undef if $feature_start >= $scaf_part_start and $feature_start <= $scaf_part_end and $feature_end > $scaf_part_end;
      $ctg_index--;
   }
   my ($ctg_start, $ctg_end) = split /-/, $ctg_intervals[$ctg_index], 2;
   
   #Convert the coordinates:
   my ($feature_ctg_start, $feature_ctg_end);
   if ($ctg_orientations[$ctg_index] eq "-") {
      $feature_ctg_start = $ctg_end - $feature_end + $scaf_part_start;
      $feature_ctg_end = $ctg_end - $feature_start + $scaf_part_start;
      $strand = $strand eq "-" ? "+" : "-";
   } else {
      $feature_ctg_start = $feature_start + $ctg_start - $scaf_part_start;
      $feature_ctg_end = $feature_end + $ctg_start - $scaf_part_start;
   }
   my $ctg = $ctg_list[$ctg_index];
   print STDERR join("\t", $feature_start, $feature_end, $feature_ctg_start, $feature_ctg_end, $ctg_start, $ctg_end, $scaf_part_start, $scaf_part_end, $ctg_orientations[$ctg_index], $ctg_index, $ctg), "\n" if $debug > 2;
   return "${ctg}:${feature_ctg_start}-${feature_ctg_end}(${strand})";
}

#Set up the hashes for translating between coordinate spaces:
my %ctg_intervals = ();
my %scaf_intervals = ();
my %scaf_ctg_map = ();
my %ctg_scaf_map = ();
my %ctg_orientations = ();
my %revised_gff3_headers = (); #An extra hash to store the revised GFF3 headers
my %scaffolded_ctgs = (); #An extra hash for quick identification of feed-through sequences
my %broken_features = (); #A hash to keep track of broken features (and check for children)

#Read through the AGP and fill the aforementioned hashes:
print STDERR "Reading AGP\n";
my $agp_fh;
open($agp_fh, "<", $input_agp);
while (!eof($agp_fh)) {
   my $line = <$agp_fh>;
   next if $line =~ /^#/;
   chomp $line;
   my @agp_parts = split /\t/, $line;
   next if $agp_parts[4] eq "N" or $agp_parts[4] eq "U"; #Skip gap lines
   
   #Collect the important information for superscaffolds:
   my $scaf = $agp_parts[0];
   my $scaf_start = $agp_parts[1];
   my $scaf_end = $agp_parts[2];
   my $ctg = $agp_parts[5];
   my $ctg_start = $agp_parts[6];
   my $ctg_end = $agp_parts[7];
   my $ctg_orientation = $agp_parts[8];
   my $ctg_interval = join("-", $ctg_start, $ctg_end);
   my $scaf_interval = join("-", $scaf_start, $scaf_end);
   
   #Fill the hashes:
   print STDERR "Multiple segments contain the same contig ${ctg} for scaffold ${scaf}\n" if exists($scaffolded_ctgs{$ctg}) and $debug;
   $scaffolded_ctgs{$ctg} = 1;
   $scaf_ctg_map{$scaf} = [] unless exists($scaf_ctg_map{$scaf});
   push @{$scaf_ctg_map{$scaf}}, $ctg;
   $ctg_scaf_map{$ctg} = $scaf; #Assumes that a contig only ever maps to one scaffold
   $ctg_intervals{$scaf} = [] unless exists($ctg_intervals{$scaf});
   push @{$ctg_intervals{$scaf}}, $ctg_interval;
   $scaf_intervals{$scaf} = [] unless exists($scaf_intervals{$scaf});
   push @{$scaf_intervals{$scaf}}, $scaf_interval;
   $ctg_orientations{$scaf} = [] unless exists($ctg_orientations{$scaf});
   push @{$ctg_orientations{$scaf}}, $ctg_orientation;
}
close($agp_fh);
print STDERR "Done reading AGP, found ", scalar(keys %scaf_ctg_map), " scaffolds\n";

#Quickly establish the revised GFF3 headers:
if ($invert) { #Make potentially truncated unscaffolded headers
   for my $scaf (keys %scaf_ctg_map) {
      my $num_scaffolded_ctgs = scalar(@{$scaf_ctg_map{$scaf}});
      for (my $i = 0; $i < $num_scaffolded_ctgs; $i++) {
         my $ctg = $scaf_ctg_map{$scaf}[$i];
         my ($ctg_start, $ctg_end) = split /-/, $ctg_intervals{$scaf}[$i], 2;
         $revised_gff3_headers{$ctg} = join(" ", "##sequence-region", $ctg, $ctg_start, $ctg_end);
      }
   }
} else { #Make scaffolded headers
   for my $scaf (keys %scaf_ctg_map) {
      my $num_scaffolded_ctgs = scalar(@{$scaf_intervals{$scaf}});
      my ($scaf_start, $ctg1_end) = split /-/, $scaf_intervals{$scaf}[0], 2;
      my ($lastctg_start, $scaf_end) = split /-/, $scaf_intervals{$scaf}[$num_scaffolded_ctgs-1], 2;
      $revised_gff3_headers{$scaf} = join(" ", "##sequence-region", $scaf, $scaf_start, $scaf_end);
   }
}

#Read through the GFF3, and modify the chromosome/scaffold, start, end, and
# orientation fields according to the hashes derived from the AGP:
my @gff_headers = ();
print STDERR "Converting input GFF3 to new coordinate space\n" unless $input_bed ne "";

#Allow for BED input and output:
my $broken_ext;
$broken_ext = "gff3?" if $input_gff ne "";
$broken_ext = "bed" if $input_bed ne "";

#Open the output GFF3 for broken features if one was provided:
my $broken_fh;
if ($broken_gff3 =~ /\.${broken_ext}$/) {
   open($broken_fh, ">", $broken_gff3);
}

my $headers_needed = 0; #Only output headers once
my $gff_fh;
open($gff_fh, "<", $input_gff) if $input_gff ne "";
open($gff_fh, "<", $input_bed) if $input_bed ne "" and $input_gff eq "";
while (!eof($gff_fh)) {
   my $line = <$gff_fh>;
   chomp $line;
   #Accumulate headers
   if ($line =~ /^##/) {
      print $broken_fh $line, "\n" if $broken_gff3 =~ /\.${broken_ext}$/ and $input_gff ne "";
      push @gff_headers, $line unless $line =~ /sequence-region/;
      if ($line =~ /sequence-region/) {
         my ($headertype, $seqid, $seq_start, $seq_end) = split /\s+/, $line, 4;
         #Keep the header if it's missing from the AGP:
         push @gff_headers, $line unless (!$invert and exists($scaffolded_ctgs{$seqid})) or ($invert and exists($scaf_ctg_map{$seqid}));
      }
      next;
   }
   #Add the revised sequence-region headers once we're onto the records:
   if ($line !~ /^##/ and $headers_needed == 0) {
      if (scalar(@gff_headers) > 0) {
         #Output the accumulated headers:
         print STDOUT join("\n", @gff_headers), "\n" if $input_gff ne "";
#         print $broken_fh join("\n", @gff_headers), "\n" if $broken_gff3 =~ /\.gff3?$/;
      }
      #Output the sorted revised headers:
      for my $scaf (sort keys %revised_gff3_headers) {
         print STDOUT $revised_gff3_headers{$scaf}, "\n" if $input_gff ne "";
#         print $broken_fh $revised_gff3_headers{$scaf}, "\n" if $broken_gff3 =~ /\.gff3?$/;
      }
      @gff_headers = ();
      $headers_needed++;
   }
   #Parse important elements out of the GFF3 records:
   my @gff_parts = split /\t/, $line;
   my $seqid = $gff_parts[0];
   #Feed-through when not found in AGP:
   unless ((!$invert and exists($scaffolded_ctgs{$seqid})) or ($invert and exists($scaf_ctg_map{$seqid}))) {
      print STDOUT $line;
      print STDOUT "#${seqid} not mapped to a chromosome\n" if $debug;
      next;
   }
   #Parse a record:
   my ($feature_start, $feature_end, $feature_orientation, $id, $parent);
   ($id, $parent) = ('', '');
   if ($input_gff ne "") { #Parse GFF record
      $feature_start = $gff_parts[3];
      $feature_end = $gff_parts[4];
      $feature_orientation = $gff_parts[6];
#      chomp $gff_parts[8];
      my @feature_tags = split /;/, $gff_parts[8];
      for my $tag (@feature_tags) {
         my ($tag_name, $tag_value) = split /=/, $tag, 2;
         $id = $tag_value if $tag_name eq "ID";
         $parent = $tag_value if $tag_name eq "Parent";
      }
   } elsif ($input_bed ne "") { #Parse BED record
      $feature_start = $gff_parts[1]+1;
      $feature_end = $gff_parts[2];
      $feature_orientation = scalar(@gff_parts) >= 6 ? $gff_parts[5] : "+";
#      chomp $gff_parts[$#gff_parts];
   }
   
   #Perform the conversions:
   if ($invert) { #Scaffold to contig
      my $scaf_coordinates = "${feature_start}-${feature_end}(${feature_orientation})";
      #Perhaps check for existence of this scaffold in the keys first?
      my $ctg_coordinates = scaffoldToContig($scaf_coordinates, $scaf_ctg_map{$seqid}, $ctg_intervals{$seqid}, $scaf_intervals{$seqid}, $ctg_orientations{$seqid}, $debug);
      #Check if this feature is broken, if so output to broken GFF3:
      unless (defined($ctg_coordinates) and !exists($broken_features{$parent})) {
         print $broken_fh $line, "\n" if $broken_gff3 =~ /\.${broken_ext}$/;
         print STDERR join("\t", $seqid, $scaf_coordinates, $id, defined($ctg_coordinates) ? "Valid coords" : "Invalid coords", exists($broken_features{$parent}) ? $parent : "Does not have broken parent"), "\n" if $debug > 1;
         $broken_features{$id} = 1 unless $id eq "";
         next; #Skip the print outside the if when it's broken
      }
      #Else, if normal, break the coordinates apart and substitute them in:
      my ($ctg, $ctg_interval_strand) = split /:/, $ctg_coordinates, 2;
      my ($ctg_interval, $strand) = split /\(/, $ctg_interval_strand, 2;
      $strand = substr($strand, 0, -1);
      my ($ctg_start, $ctg_end) = split /-/, $ctg_interval, 2;
      $gff_parts[0] = $ctg;
      if ($input_gff ne "") {
         $gff_parts[3] = $ctg_start;
         $gff_parts[4] = $ctg_end;
         $gff_parts[6] = $strand;
      } elsif ($input_bed ne "") {
         $gff_parts[1] = $ctg_start-1;
         $gff_parts[2] = $ctg_end;
         $gff_parts[5] = $strand if scalar(@gff_parts) >= 6;
      }
   } else { #Contig to scaffold
      my $ctg_coordinates = "${seqid}:${feature_start}-${feature_end}(${feature_orientation})";
      my $scaf = $ctg_scaf_map{$seqid};
      my $scaf_coordinates = contigToScaffold($ctg_coordinates, $scaf_ctg_map{$scaf}, $ctg_intervals{$scaf}, $scaf_intervals{$scaf}, $ctg_orientations{$scaf}, $debug);
      #Check if this feature is broken, if so output to broken GFF3:
      unless (defined($scaf_coordinates) and !exists($broken_features{$parent})) {
         print $broken_fh $line, "\n" if $broken_gff3 =~ /\.${broken_ext}$/;
         print STDERR join("\t", $seqid, $ctg_coordinates, $id, defined($scaf_coordinates) ? "Valid coords" : "Invalid coords", exists($broken_features{$parent}) ? $parent : "Does not have broken parent"), "\n" if $debug > 1;
         $broken_features{$id} = 1 unless $id eq "";
         next; #Skip the print outside the if when it's broken
      }
      my ($scaf_interval, $strand) = split /\(/, $scaf_coordinates, 2;
      $strand = substr($strand, 0, -1);
      my ($scaf_start, $scaf_end) = split /-/, $scaf_interval, 2;
      $gff_parts[0] = $scaf;
      if ($input_gff ne "") {
         $gff_parts[3] = $scaf_start;
         $gff_parts[4] = $scaf_end;
         $gff_parts[6] = $strand;
      } elsif ($input_bed ne "") {
         $gff_parts[1] = $scaf_start-1;
         $gff_parts[2] = $scaf_end;
         $gff_parts[5] = $strand if scalar(@gff_parts) >= 6;
      }
   }

   #Print the modified GFF line:
   print STDOUT join("\t", @gff_parts), "\n";
}
close($gff_fh);
close($broken_fh) if $broken_gff3 =~ /\.${broken_ext}$/;
print STDERR "Done converting GFF3 to new coordinate space\n";

exit 0;
