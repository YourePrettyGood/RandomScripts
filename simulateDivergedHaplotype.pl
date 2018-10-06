#!/usr/bin/env perl

use warnings;
use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

############################################################################################
# simulateDivergedHaplotype.pl                                                             #
# Usage:                                                                                   #
#  simulateDivergedHaplotype.pl [-i reference FASTA] [-o output FASTA] [-n] [-g geom param]#
#   <% divergence>                                                                         #
#                                                                                          #
# Arguments:                                                                               #
#  -i,--input_haplotype         FASTA of the first haplotype to use as a reference for     #
#                               homozygous sites (optional, defaults to STDIN)             #
#  -o,--output_haplotype        Output FASTA of the new diverged haplotype to a file       #
#                               (optional, default outputs to STDOUT)                      #
#  -n,--indels                  If present, this option adds indels at 4% of the divergence#
#                               as per Drosophila estimates)                               #
#                               (flag, only include if you do want indels)                 #
#  -g,--indel_geom              Parameter value for the geometric distribution used to     #
#                               model indel length distribution                            #
#                               (optional, decimal between 0 and 1, default: 0.1)          #
#  % divergence                 Required parameter specifying the percent divergence       #
#                               desired between the reference/input haplotype and the      #
#                               new/output diverged haplotype                              #
#                               (required, specify as decimal between 0 and 100)           #
#                                                                                          #
# Description:                                                                             #
#  simulateDivergedHaplotype.pl simulates a diverged haplotype given a reference haplotype #
#  and an amount of divergence based on a simple statistical model:                        #
#  scaffold length * % divergence / 100 = # SNPs on the scaffold                           #
#  # SNPs on scaffold * 0.04 * 1 / 2 = # deletions on the scaffold (if indels allowed)     #
#  # SNPs on scaffold * 0.04 * 1 / 2 = # insertions on the scaffold (if indels allowed)    #
#  SNPs and indels are distributed along scaffolds under a discrete uniform distribution   #
#  Indel length is distributed under a geometric distribution with parameter 0.1 by default#
#  This value was very qualitatively estimated from Rimmer et al. 2014 Nat. Gen. for humans#
############################################################################################

my $SCRIPTNAME = "simulateDivergedHaplotype.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

simulateDivergedHaplotype.pl - Simulate a haplotype diverged from a reference haplotype

=head1 SYNOPSIS

simulateDivergedHaplotype.pl [options] <% divergence>

 Options:
  --help,-h,-?             Display this help documentation
  --input_haplotype,-i     FASTA of the reference haplotype
                           (optional, default: STDIN)
  --output_haplotype,-o    Output the new diverged haplotype to this FASTA
                           file
                           (optional, default: STDOUT)
  --indels,-n              If present, this option adds indels at 4% of the
                           divergence as per Drosophila estimates, composed
                           of 5 deletions per insertion (again, per
                           Drosophila estimates)
                           (flag, only include if you DO want indels)
  --indel_geom,-g          Parameter value for the geometric distribution
                           used to model indel length distribution
                           (optional, decimal between 0 and 1, default: 0.1)
  --version,-v             Output version string

 Mandatory:
  % divergence             Required parameter specifying the percent
                           divergence desired between the reference
                           haplotype and the output haplotype
                           (required, specify as a real number between 0 and
                           100)

=head1 DESCRIPTION

simulateDivergedHaplotype.pl simulates SNPs and indels at a given rate onto a supplied
reference haplotype, thereby generating a haplotype with a specific degree of divergence
to the reference.  The length of the scaffold times the percent divergence gives the
number of SNPs introduced, the locations of which are drawn from a discrete uniform
distribution over the length of the scaffold.  The substituted bases are drawn with equal
probability from each of the alternate bases (e.g. if ref is A, Pr(C)=Pr(G)=Pr(T)=1/3).
Indels are simulated at 4% of the SNP rate.  Indel lengths are drawn from a geometric
distribution with adjustable parameter.  The default parameter value, p=0.1, was
qualitatively estimated from Rimmer et al. 2014 Nature Genetics human data.  Two logs are
generated to keep track of what and where SNPs and indels were generated.

This script assumes that the input haplotype FASTA is NOT line-wrapped.

=cut

# sampleWithoutReplacement($max, $number_of_samples, $reference_to_hash)
# Sets $number_of_samples elements of the hash passed in
# based on sampling without replacement from a discrete 
# uniform distribution from 0 to $max.
# TODO: Improve this with a bitvector instead of a hash
sub sampleWithoutReplacement($$) {
   my $max = shift @_;
   my $number_of_samples = shift @_;
   my %samples = ();
   #Based on http://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
   # which is based on Algorithm 3.4.2S of Knuth's Seminumeric Algorithms
   my ($t, $samples_drawn) = (0, 0);
   while ($samples_drawn < $number_of_samples) {
      my $u = rand();
      if (($max - $t)*$u < $number_of_samples - $samples_drawn) {
         $samples{$t} = 1;
         $samples_drawn++;
      }
      $t++;
   }
   return \%samples;
}

# insertMutations($divergence, $scaffold, $sorted_SNP_array, $sorted_indel_array)
# Returns the mutated scaffold
# Indel length simulation is done in this function
sub insertMutations($$$$$) {
   my $divergence = shift @_;
   my $scaffold = shift @_;
   my %SNPs = %{shift @_};
   my %indels = %{shift @_};
   my $indel_geometric_parameter = shift @_;

   my @int_to_nuc = ("A","C","G","T");
   
   my $mutated_scaffold = "";
   my @SNP_log = ();
   my @indel_log = ();
   my ($SNP_index, $indel_index) = (0, 0);
   my ($SNP_hash_size, $indel_hash_size) = (scalar(keys %SNPs), scalar(keys %indels));
   my $scaffold_length = length($scaffold);
   for (my $i = 0; $i < $scaffold_length; $i++) {
      my $ref = substr $scaffold, $i, 1;
      if ($SNP_index < $SNP_hash_size and defined($SNPs{$i})) { #SNP site, so mutate it
         #Note: If SNP and indel positions overlap, we mutate the SNP, but do not introduce the indel
         #Hence, the effective indel rate is slightly lower than expected.
         my @alt_nucs = grep { $_ ne uc $ref } @int_to_nuc; #Generate an array of non-ref bases
         my $new_base = $alt_nucs[int(rand(3))]; #Choose a random non-ref base
         $mutated_scaffold .= $new_base; #Mutate the ref
         push @SNP_log, join("\t", $i+1, $ref, $new_base);
         $SNP_index++;
      } elsif ($indel_index < $indel_hash_size and defined($indels{$i})) { #indel site
         #Inverse CDF method for generating random number from geometric distribution:
         my $indel_length = int(log(1-rand())/log(1-$indel_geometric_parameter));
         #Make sure we don't create a 0-length indel:
         $indel_length = int(log(1-rand())/log(1-$indel_geometric_parameter)) while $indel_length == 0;
         if (int(rand(2)) == 1) { #Choice of 1 is arbitrary, but this ensures 1/2 chance of insertion
            my $indel = "";
            for (my $j = 1; $j <= $indel_length; $j++) { #Compose the indel sequence (random nucleotides, no GC bias)
               $indel .= $int_to_nuc[int(rand(4))];
            }
            $mutated_scaffold .= substr($scaffold, $i, 1) . $indel;
            push @indel_log, join("\t", $i+1, "ins", $indel_length, $indel);
         } else { #Deletion skips $indel_length bases in ref including the current base
            if ($i+$indel_length >= $scaffold_length) { #Unfortunately, this deletion will be smaller than expected:
               push @indel_log, join("\t", $i+1, "del", $indel_length, substr($scaffold, $i));
            } else {
               push @indel_log, join("\t", $i+1, "del", $indel_length, substr($scaffold, $i, $indel_length));
            }
            $i += $indel_length-1; #The loop post-increment takes care of the remaining 1
         }
         $indel_index++;
      } else { #Unchanged site
         $mutated_scaffold .= substr $scaffold, $i, 1;
      }
   }
   return ($mutated_scaffold, \@SNP_log, \@indel_log);
}

#Parse options:
my $help = 0;
my $man = 0;
my $ref_path = '';
my $out_path = '';
my $indels = 0;
my $indel_geometric_parameter = 0.1;
my $dispversion = 0;
GetOptions('input_haplotype|i=s' => \$ref_path, 'output_haplotype|o=s' => \$out_path, 'indels|n' => \$indels, 'indel_geom|g=f' => \$indel_geometric_parameter, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $percent_divergence = 1.0;
if (scalar@ARGV < 1) { #Missing percent divergence parameter
   print STDERR "Missing the percent divergence parameter, exiting.\n";
   exit 2;
} else {
   $percent_divergence = $ARGV[0];
}
my $divergence = $percent_divergence/100;

#Open input and output files:
my ($ref, $out);
if ($ref_path =~ /\.gz$/) {
   $ref = new IO::Uncompress::Gunzip $ref_path or die "Failed to open Gzipped reference FASTA file $ref_path due to error $GunzipError\n";
} elsif ($ref_path eq '') {
   open $ref, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for reference FASTA due to error $!\n";
} else {
   open $ref, "<", $ref_path or die "Failed to open reference FASTA file due to error $!\n";
}
if ($out_path =~ /\.gz$/) {
   $out = new IO::Compress::Gzip $out_path or die "Failed to create Gzipped output FASTA file $out_path due to error $GzipError\n";
} elsif ($out_path eq '') {
   open $out, ">&", \*STDOUT or die "Failed to duplicate STDOUT file handle for output FASTA due to error $!\n";
} else {
   open $out, ">", $out_path or die "Failed to open output FASTA file due to error $!\n";
}
#Generate output SNP and indel log files:
$out_path =~ s/\.\w+(.gz)?$/_/;
my $log_prefix = $out_path;
my $SNP_log_path = $log_prefix . "SNPs.log";
my $indel_log_path = $log_prefix . "indels.log";
open SNPLOG, ">", $SNP_log_path or die "Failed to open output SNP log file $SNP_log_path due to error $!\n";
open INDELLOG, ">", $indel_log_path or die "Failed to open output indel log file $indel_log_path due to error $!\n";

#SNP rate relative to indel rate:
my $indel_rate_fold_lower = 25;
#Iterate over the scaffolds:
#Assumes FASTA is not line-wrapped
my $scaffold_name = ""; #Keep track of the scaffold name
while (my $line = <$ref>) {
   if ($line =~ />/) {
      print $out $line;
      $scaffold_name = $line;
      chomp $scaffold_name;
      $scaffold_name =~ s/>//;
      next;
   } else {
      chomp $line;
      my %SNPs = ();
      my %indels = ();
      my $scaffold_length = length($line);
      my $num_SNPs = int($scaffold_length*$divergence);
      #Simulate SNP and indel positions from a random discrete uniform distribution:
      %SNPs = %{sampleWithoutReplacement($scaffold_length-1, $num_SNPs)};
      %indels = %{sampleWithoutReplacement($scaffold_length-1, int($num_SNPs/$indel_rate_fold_lower))} if $indels != 0;
      #Diagnostic:
      print STDERR "Expecting to output ", scalar(keys %SNPs), " SNPs and ", scalar(keys %indels), " indels.\n";
      #End diagnostic
      #Output the desired mutated scaffold:
      my ($mutated_scaffold, $SNP_log_reference, $indel_log_reference) = insertMutations($divergence, $line, \%SNPs, \%indels, $indel_geometric_parameter);
      print $out $mutated_scaffold, "\n";
      #Output the SNP and indel logs:
      for my $SNP (@{$SNP_log_reference}) {
         print SNPLOG $scaffold_name, "\t", $SNP, "\n";
      }
      for my $indel (@{$indel_log_reference}) {
         print INDELLOG $scaffold_name, "\t", $indel, "\n";
      }
      #Diagnostic:
      print STDERR "Really output ", scalar@{$SNP_log_reference}, " SNPs and ", scalar@{$indel_log_reference}, " indels.\n";
      #End diagnostic
   }
}

#Close input and output files:
close $ref;
close $out;
close SNPLOG;
close INDELLOG;

exit 0;
