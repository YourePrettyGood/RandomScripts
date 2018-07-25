#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#Create a GFF3 file of minimal gene annotations given the output
# of exonerate mapping proteins to a genome.
#The "minimal" gene annotation does not include any UTRs, so each
# exon record is identical to its CDS record, and each gene record
# is similar if not identical to its mRNA record(s).

=pod

=head1 NAME

createGFFfromExonerate.pl - Create a minimal GFF3 from exonerate output (must have VULGAR line)

=head1 SYNOPSIS

createGFFfromExonerate.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --vulgar_path,-v      Path to exonerate output files
                        (must contain VULGAR alignment strings)
  --wonky_gff,-w        Path to output GFF3 of wonky genes
                        (optional, recommended to output)
  --genome_fai,-f       Path to genome FASTA index (.fai) file
                        (required, used for GFF header lines)
  --prefix,-p           Prefix for genes and transcripts if no gene or
                        transcript IDs found in Query: line of exonerate
                        output

=head1 DESCRIPTION

This script creates a minimal GFF3 gene annotation for a genome based
on exonerate alignments of proteins to the genome.  We use VULGAR
alignment strings to define the exons and CDSes, and use the protein's
FASTA header to define the gene and mRNA IDs corresponding to the
query protein (e.g. the FBgn and FBtr corresponding to an FBpp).

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $vulgar_path = "";
my $wonky_path = "";
my $fai_path = "";
my $prefix = "my";
GetOptions('vulgar_path|v=s' => \$vulgar_path, 'wonky_gff|w=s' => \$wonky_path, 'genome_fai|f=s' => \$fai_path, 'prefix|p=s' => \$prefix, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

#Open the genome FASTA index (.fai) file and get set for the GFF3
# ##sequence-region header lines:
die "Missing .fai file, cannot proceed (or will produce invalid GFF3 files)" if $fai_path eq "";
open(FAI, "<", $fai_path) or die "Cannot open .fai file, cannot proceed";
my %scaf_headers = ();
while (my $fai_line = <FAI>) {
   chomp $fai_line;
   my @fai_arr = split /\t/, $fai_line;
   $scaf_headers{$fai_arr[0]} = join(" ", "##sequence-region", $fai_arr[0], "1", $fai_arr[1]);
}
close(FAI);

#Open the exonerate alignment VULGAR file:
unless(open(VULGAR, "<", $vulgar_path)) {
   print STDERR "Error opening exonerate alignment VULGAR file.\n";
   exit 3;
}

sub findMaxRange($$) {
   my @transcripts = @{shift @_};
   my %transcript_ranges = %{shift @_};
   my $min = 100000000;
   my $max = 0;
   my $scaf = '';
   my $strand = '';
   my $wonky = 0;
   for my $transcript (@transcripts) {
      my @range = @{$transcript_ranges{$transcript}};
      if ($range[0] ne $scaf and $scaf ne '') {
         print STDERR "Gene with transcript ${transcript} has mRNAs on multiple scaffolds\n";
         $wonky = 1;
         next; #Do not include second scaffold hit's positions in min and max
      } elsif ($range[1] ne $strand and $strand ne '') {
         print STDERR "Gene with transcript ${transcript} has mRNAs on multiple strands\n";
         $wonky = 1;
         next; #Do not include opposite strand positions in min and max
      }
      $scaf = $range[0];
      $strand = $range[1];
      $min = $range[2] unless $range[2] > $min;
      $max = $range[3] unless $range[3] < $max;
   }
   return([$wonky, $scaf, $strand, $min, $max]);
}

sub vulgarToCDSes($$$) {
   my @vulgaroptuples = @{shift @_};
   my @intervals = ();
   my $num_optuples = scalar@vulgaroptuples/3;
   my $interval_start = shift @_;
   my $strand = shift @_;
   my @phases = ();
   my $phase = 0;
   push @phases, $phase; #First exon's phase is always 0
   my $phase_shift = 0; #Shift the phase when frameshifts occur
   my $exonstart = 0; #Flag to trigger saving of phase upon starting a new exon
   my @CDSes = ();
   my $interval_end = $interval_start;
   for (my $i = 0; $i < $num_optuples; $i++) {
      my $op = $vulgaroptuples[3*$i];
      my $prot_len = $vulgaroptuples[3*$i+1];
      my $dna_len = $vulgaroptuples[3*$i+2];
      if ($op eq "M") {
         if ($exonstart) {
            push @phases, (($phase + $phase_shift) % 3);
            $exonstart = 0;
            $phase = 0;
         }
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
      } elsif ($op eq "G") {
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
         #We don't need to do anything for deletions, as $dna_len is 0 for them
      } elsif ($op eq "S") {
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
         $phase = $dna_len;
      } elsif ($op eq "F") {
         #Do we ignore frameshifts and output the interval without it?
         #Or do we report incorrect phase while including frameshifted bases?
         #We use a phase shift to report the correct phase, and output
         #intervals including the frameshifting bases
         $phase_shift += $dna_len;
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
      } elsif ($op eq "5") {
         push @intervals, join("\t", $interval_start, $interval_end-1) if $strand eq "+";
         push @intervals, join("\t", $interval_end+1, $interval_start) if $strand eq "-";
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
         $interval_start = $interval_end;
      } elsif ($op eq "I" or $op eq "3") {
         $interval_end += $dna_len if $strand eq "+";
         $interval_end -= $dna_len if $strand eq "-";
         $interval_start = $interval_end;
         $exonstart = 1 if $op eq "3";
      } else {
         die "Unknown VULGAR op ${op}\n";
      }
   }
   push @intervals, join("\t", $interval_start, $interval_end-1) if $vulgaroptuples[-3] eq "M" and $strand eq "+";
   push @intervals, join("\t", $interval_end+1, $interval_start) if $vulgaroptuples[-3] eq "M" and $strand eq "-";
   
   for (my $i = 0; $i <= $#intervals; $i++) {
      #Columns 4-8 of GFF line:
      die "Strand missing" unless defined($strand);
      die "Phase missing for CDS ${i}" unless $i <= $#phases;
      push @CDSes, join("\t", $intervals[$i], ".", $strand, $phases[$i]);
   }
   return(\@CDSes);
}

my %CDS_intervals = ();
my %CDS_scaffolds = ();
my %prot_to_gene = (); #Maps FBpps to FBgns, should be onto, maybe not 1:1
my %prot_to_transcript = (); #Maps FBpps to FBtrs, should be 1:1 and onto
my %transcripts = (); #Maps FBgns to FBtrs, should be onto, likely not 1:1
my %exons = (); #Maps FBtrs to lists of GFF strings for exons and CDSes
my %transcript_strings = (); #Maps FBtrs to GFF strings for mRNAs
my %transcript_ranges = (); #Maps FBtrs to start and end coordinates (for later sorting and determining start and end of genes)
my %wonky_transcripts = (); #For abnormal cases, e.g. multiple equal best alignments
my %genes = (); #Maps FBgns to GFF strings for genes
my %wonky_genes = (); #For abnormal cases, e.g. gene products map to multiple chromosomes
my %gene_ranges = (); #Maps FBgns to start and end coordinates for sorting
my %wonky_gene_ranges = (); #For abnormal cases
while (my $line = <VULGAR>) {
   chomp $line;
   next unless $line =~ /(Query|vulgar):/;
   if ($line =~ /Query:\s+(\S+)\s+/) {
      my $protein_id = $1;
      my $gene_id;
      my $transcript_id;
      if ($line =~ /parent=(\S+);/) { #This is sort of hacky, but works for FlyBase protein headers
         my $parents = $1;
         my @parentarr = split /,/, $parents;
         #The following line is specific to FlyBase proteins
         #die "More than two parents for protein ${protein_id}: $parents" unless $parentarr[0] =~ /^FBgn/ and $parentarr[1] =~ /^FBtr/;
         die "More than two parents for protein ${protein_id}: $parents" if scalar(@parentarr) > 2; #Protein can have gene and transcript parents, so > 2 parents means > 1 gene
         $gene_id = $parentarr[0];
         $transcript_id = $parentarr[1];
      } else { #Just force them to have identical ID numbers, but extra prefixes
         $gene_id = "${prefix}gn${protein_id}";
         $transcript_id = "${prefix}tr${protein_id}";
      }
      $prot_to_transcript{$protein_id} = $transcript_id;
      $prot_to_gene{$protein_id} = $gene_id;
      $transcripts{$gene_id} = [] unless exists($transcripts{$gene_id});
      push @{$transcripts{$gene_id}}, $transcript_id;
   } else {
      my @vulgararr = split /\s+/, $line;
      shift @vulgararr; #"vulgar:"
      my $prot_id = shift @vulgararr;
      shift @vulgararr; #Residue start
      shift @vulgararr; #Residue end
      shift @vulgararr; #Meaningless protein strand
      my $prot_scaf = shift @vulgararr;
      my $prot_start = shift @vulgararr;
      my $prot_end = shift @vulgararr;
      my $prot_strand = shift @vulgararr;
      $prot_start += 1 if $prot_strand eq "+";
      shift @vulgararr; #Alignment score
      #Flag this mRNA as wonky if it already exists:
      $wonky_transcripts{$prot_to_transcript{$prot_id}} = 1 if exists($transcript_strings{$prot_to_transcript{$prot_id}});
      #Construct mRNA GFF string:
      $transcript_strings{$prot_to_transcript{$prot_id}} = join("\t", $prot_scaf, "Exonerate", "mRNA", $prot_start, $prot_end+3, ".", $prot_strand, ".", "ID=" . $prot_to_transcript{$prot_id} . ";Parent=" . $prot_to_gene{$prot_id}) if $prot_strand eq "+";
      $transcript_strings{$prot_to_transcript{$prot_id}} = join("\t", $prot_scaf, "Exonerate", "mRNA", $prot_end-3, $prot_start, ".", $prot_strand, ".", "ID=" . $prot_to_transcript{$prot_id} . ";Parent=" . $prot_to_gene{$prot_id}) if $prot_strand eq "-";
      $transcript_ranges{$prot_to_transcript{$prot_id}} = [$prot_scaf, $prot_strand, $prot_start, $prot_end+3] if $prot_strand eq "+";
      $transcript_ranges{$prot_to_transcript{$prot_id}} = [$prot_scaf, $prot_strand, $prot_end-3, $prot_start] if $prot_strand eq "-";
      #Identify CDSes from VULGAR string:
      my @CDSes = @{vulgarToCDSes(\@vulgararr, $prot_start, $prot_strand)};
      #Construct exon and CDS GFF strings:
      $exons{$prot_to_transcript{$prot_id}} = [] unless exists($exons{$prot_to_transcript{$prot_id}});
      for my $CDS (@CDSes) {
         my $CDS_string = join("\t", $prot_scaf, "Exonerate", "CDS", $CDS, "Parent=" . $prot_to_transcript{$prot_id});
         my $exon_string = join("\t", $prot_scaf, "Exonerate", "exon", $CDS, "Parent=" . $prot_to_transcript{$prot_id});
         push @{$exons{$prot_to_transcript{$prot_id}}}, join("\n", $exon_string, $CDS_string);
      }
   }
}

close(VULGAR);

#Generate gene GFF strings with ranges based on the min start and max end of
# component transcripts:
for my $fbgn (keys %transcripts) {
   my $wonky = 0;
   for my $fbtr (@{$transcripts{$fbgn}}) {
      if (exists($wonky_transcripts{$fbtr})) {
         $wonky = 1;
         last;
      }
   }
   my @gene_range = @{findMaxRange(\@{$transcripts{$fbgn}}, \%transcript_ranges)};
   #First element is a boolean indicating whether the gene is wonky or not
   my $gene_wonky = shift @gene_range;
   my $gene_string = join("\t", $gene_range[0], "Exonerate", "gene", $gene_range[2], $gene_range[3], ".", $gene_range[1], ".", "ID=${fbgn}");
   $genes{$fbgn} = $gene_string unless $wonky or $gene_wonky;
   $wonky_genes{$fbgn} = $gene_string if $wonky or $gene_wonky;
   $gene_ranges{$fbgn} = \@gene_range unless $wonky or $gene_wonky;
   $wonky_gene_ranges{$fbgn} = \@gene_range if $wonky or $gene_wonky;
}

#Finally, sort the genes by scaffold and position, and output the GFF:
my %genes_by_scaf = ();
#Assign genes to scaffolds:
for my $fbgn (keys %gene_ranges) {
   my @gene_range = $gene_ranges{$fbgn};
   $genes_by_scaf{$gene_range[0]} = () unless exists($genes_by_scaf{$gene_range[0]});
   push @{$genes_by_scaf{$gene_range[0]}}, $fbgn;
}
#Sort genes by start position in scaffolds:
my %gene_order_by_scaf = ();
for my $scaf (keys %genes_by_scaf) {
   @{$gene_order_by_scaf{$scaf}} = sort {$gene_ranges{$a}[2] <=> $gene_ranges{$b}[2] } @{$genes_by_scaf{$scaf}};
}

#Output the GFF header:
print STDOUT "##gff-version 3\n";
#Print ##sequence-region headers:
for my $scaf (sort keys %scaf_headers) {
   print STDOUT $scaf_headers{$scaf}, "\n";
}

#Output the GFF records:
for my $scaf (keys %gene_order_by_scaf) {
   for my $fbgn (@{$gene_order_by_scaf{$scaf}}) {
      print STDOUT $genes{$fbgn}, "\n";
      for my $fbtr (@{$transcripts{$fbgn}}) {
         print STDOUT $transcript_strings{$fbtr}, "\n";
         for my $exon (@{$exons{$fbtr}}) {
            print STDOUT $exon, "\n";
         }
      }
   }
}

if ($wonky_path ne '') {
   die "Unable to open output GFF3 for wonky genes: $wonky_path" unless open(WONKY, ">", $wonky_path);
   #Finally, sort the genes by scaffold and position, and output the GFF:
   my %wonky_genes_by_scaf = ();
   #Assign genes to scaffolds:
   for my $fbgn (keys %wonky_gene_ranges) {
      my @wonky_gene_range = $wonky_gene_ranges{$fbgn};
      $wonky_genes_by_scaf{$wonky_gene_range[0]} = () unless exists($wonky_genes_by_scaf{$wonky_gene_range[0]});
      push @{$wonky_genes_by_scaf{$wonky_gene_range[0]}}, $fbgn;
   }
   #Sort genes by start position in scaffolds:
   my %wonky_gene_order_by_scaf = ();
   for my $scaf (keys %wonky_genes_by_scaf) {
      @{$wonky_gene_order_by_scaf{$scaf}} = sort {$wonky_gene_ranges{$a}[2] <=> $wonky_gene_ranges{$b}[2] } @{$wonky_genes_by_scaf{$scaf}};
   }
   #Output the GFF header:
   print WONKY "##gff-version 3\n";
   #Optionally, we could use info from the FAI to construct sequence-region headers
   
   #Output the GFF records:
   for my $scaf (keys %wonky_gene_order_by_scaf) {
      for my $fbgn (@{$wonky_gene_order_by_scaf{$scaf}}) {
         print WONKY $wonky_genes{$fbgn}, "\n";
         for my $fbtr (@{$transcripts{$fbgn}}) {
            print WONKY $transcript_strings{$fbtr}, "\n";
            for my $exon (@{$exons{$fbtr}}) {
               print WONKY $exon, "\n";
            }
         }
      }
   }
   close(WONKY);
}

#Now we can iterate through the genome FASTA, and store the records
#my %scaffolds = ();
#my $scaffold_name = "";
#my $scaffold_sequence = "";
#while (my $line = <GENOME>) {
#   chomp $line;
#   #If we're at a header line and we've seen header lines before,
#   # output the sites from the previous scaffold (since we're on
#   # a new scaffold's header line):
#   if ($line =~ /^>/) {
#      $scaffolds{$scaffold_name} = $scaffold_sequence unless $scaffold_sequence eq "" or $scaffold_name eq "";
#      my $scaffold_name_line = substr $line, 1; #Get rid of the prefixed ">"
#      my @scaffold_name_parts = split /\s+/, $scaffold_name_line;
#      $scaffold_name = $scaffold_name_parts[0];
#      $scaffold_sequence = ""; #Clear out the old sequence
#   } else { #Sequence line
#      $scaffold_sequence .= $line;
#   }
#}
#$scaffolds{$scaffold_name} = $scaffold_sequence;
#close(GENOME);

#sub printCDS($$$) {
#   my $prot_id = shift @_;
#   my $sequence = shift @_;
#   my @intervals = @{shift @_;};
#   my $CDS = "";
#   for my $interval (@intervals) {
#      my ($start, $end) = split /-/, $interval, 2;
#      if ($start <= $end) { # + strand
#         my $segment = substr $sequence, $start-1, $end-$start+1;
#         $CDS .= uc($segment);
#      } else { # - strand
#         my $segment = substr $sequence, $end-1, $start-$end+1;
#         my $revcompseq = reverse uc($segment);
#         $revcompseq =~ tr/[ACGTRYSWKMBDHVN]/[TGCAYRSWMKVHDBN]/;
#         $CDS .= $revcompseq;
#      }
#   }
#   print STDOUT ">", $prot_id, "\n", $CDS, "\n";
#}

#for my $CDS (keys %CDS_scaffolds) {
#   my $scaffold = $CDS_scaffolds{$CDS};
#   my $intervals = $CDS_intervals{$CDS};
#   die "Scaffold ${scaffold} not found in input FASTA.\n" unless exists($scaffolds{$scaffold});
#   my $scaf_seq = $scaffolds{$scaffold};
#   print STDERR $CDS, " ", $scaffold, " ", join(":", @{$intervals}), "\n"; #Diagnostic
#   printCDS($CDS, $scaf_seq, $intervals);
#}

exit 0;
