# RandomScripts
Collection of arbitrary (perhaps useful) programs and scripts, and occasional one-liners

Most scripts that I've written have a `-h` option, so take a look, and let me know if anything is unclear! (Awk scripts here do not have any flags, though)

## Categories
1. [De novo assembly-related](README.md#de-novo-assembly-related-scripts)
1. [Annotation-related](README.md#annotation-related-scripts)
1. [Genome-wide statistics](README.md#genome-wide-statistics-programsscripts)
1. [MSA-related](README.md#msa-related-scripts)
1. [Illumina data parsing](README.md#illumina-data-parsing-scripts)

## Scripts that are not mine:
1. `barcode_splitter.py`

This script is needed by `divideConquerParser.sh`, as it is the main parsing workhorse.  divideConquerParser.sh just preprocesses the files to split into `n` parts, run in parallel on `n` cores, and merge the results back together.

## One-liners:

### Convert FASTQ to unwrapped FASTA (FASTQ can be gzipped, or uncompressed):

`gzip -dcf [FASTQ file] | awk 'NR%4==1||NR%4==2' | sed 's/@/>/' > [FASTA file]`

### Output lengths of longest `N` contigs or scaffolds from an assembly:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | awk '/^>/{scafname=$1;}!/^>/{print substr(scafname, 2)"\t"length($0);}' | sort -k2,2nr | head -n[N]`

**Note: FASTX toolkit is available from [Hannon Lab at CSHL](http://hannonlab.cshl.edu/fastx_toolkit/)**

### Split scaffolds into contig at N gaps of length >= 1:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | awk '/^>/{header=$0;}!/^>/{split($0, ctgs, /[N]+/); for (ctg in ctgs) {print header"_ctg"ctg; print ctgs[ctg];};}' > [unwrapped contigs FASTA]`

Just adjust the regex in the `split()` call if you want to split at gaps of specified minimum length.

### Generate a BED of gaps from a genome FASTA:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | perl -e 'use warnings; use strict; my $scafname = ""; while (my $line = <>) { chomp $line; if ($line =~ /^>/) {$scafname = substr $line, 1;} else { my $seq = $line; print join("\t", $scafname, $-[0], $+[0]+1), "\n" while $seq =~ /[^ACGTacgt]+/g;};}' > [gaps BED file]`

This outputs BED intervals for any regions with characters other than ACGT, so degenerate bases will be included in these intervals. Just adjust the regex if you want different gap detection behaviour (e.g. `/[Nn]+/g` if you only want to consider N gaps). The regex I used above was intended for use with a haploid reference assembly for the genome size estimation pipeline from Davey et al. (2016) G3 [doi:10.1534/g3.115.023655](https://dx.doi.org/10.1534/g3.115.023655)

## De novo assembly-related scripts:

### `NX.pl`

Calculates some standard contiguity statistics from an assembly FASTA.  By default, it calculates the N50 and L50 (length-weighted median contig/scaffold size and number, respectively), and all calls to it also indicate the total assembly size, number of contigs/scaffolds in the assembly, length of longest and shortest contig/scaffold, and average contig/scaffold length.  Options include changing the quantile from 50 (e.g. set quantile to 90 to calculate N90 and L90, setting it to 50 is the same as default), and specifying the genome size (e.g. in the default case, to calculate an NG50).

Multiple quantiles can be evaluated by specifying a comma-separated (no spaces) list of quantiles as input.

You can pass `-` as the FASTA file path, and the script will read the FASTA from STDIN, so it can easily be incorporated into piped one-liners.

Also, interestingly, with some upstream fiddling to get reads into FASTA format, you can calculate the NRX statistics (like NR25 and NR40) with this script, although it's fairly slow for large numbers of reads.  We're limited by Perl's sort routine's time complexity, unfortunately.

Example usage:

Calculate N50, L50, and other stats for Drosophila melanogaster FlyBase release 6.13 assembly:

`NX.pl dmel-all-chromosome-r6.13.fasta`

Calculate N90, L90, etc. for the same assembly:

`NX.pl dmel-all-chromosome-r6.13.fasta 90`

Calculate the NG50, LG50, etc., assuming G=130,000,000 bp:

`NX.pl dmel-all-chromosome-r6.13.fasta 50 130000000`

Calculate N50 on a gzipped assembly (e.g. from 10x Genomics' Supernova mkoutput):

`gzip -dc MySupernovaAsm.fa.gz | NX.pl -`

or

`NX.pl <(gzip -dc MySupernovaAsm.fa.gz)`

or

`NX.pl MySupernovaAsm.fa.gz`

though piping and command substitution seem to be faster in most cases.

Calculate the N10, N50, and N90 for the Dmel assembly used above:

`NX.pl dmel-all-chromosome-r6.13.fasta 10,50,90`

### `manualScaffold.pl`

Usage:

`manualScaffold.pl -i [path to unscaffolded FASTA] [options] <configuration string>`

This was originally developed as a quick way to manually scaffold contigs into chromosome arms.  You must know *a priori* what the order and orientation of contigs needs to be, as you specify a configuration string, or supply an AGP file to dictate how the script sews the contigs together, and how to label the resultant scaffolds.

If you wish to read the input unscaffolded FASTA from STDIN, just omit the `-i` option.

The `-u` option prints contigs that were not scaffolded into the output FASTA as well.  It just saves you from having to write out a very long configuration string or AGP file.

It has now been adapted to scaffold based on an [AGP version 2.0](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file with the `-a` option, which properly treats `N` and `U` gap records for length, although this script ignores columns 2, 3, 7, and 8 (so it takes the entirety of the contig even if the AGP says to take a substring).

The original configuration string format isn't the most intuitive, but it was easy to come up with.  Defined using the [augmented BNF grammar](https://en.wikipedia.org/wiki/Augmented_Backus%E2%80%93Naur_form):

`config-string = scaffold*(scaffold-delim scaffold)`

`scaffold-delim = "=>"`

`scaffold = scaffold-name record-delim contig*(contig-delim contig)`

`record-delim = ":"`

`contig-delim = "->"`

`contig = contig-name [revcomp]`

`revcomp = "*"`

`contig-name = *(ALPHA / DIGIT / "_" / "|")`

So for example, if I had 3 contigs generated by Quiver polishing of FinisherSC-improved contigs, and I knew that they should be given the scaffold name chr1 and joined as Segkk0|quiver's end with Segkk1|quiver's beginning, and Segkk1|quiver's end with Segkk2|quiver's end, the configuration string would be:

`chr1:Segkk0|quiver->Segkk1|quiver->Segkk2|quiver*`

If there were more scaffolds to make, the configuration string might look like:

`chr1:Segkk0|quiver->Segkk1|quiver->Segkk2|quiver*=>chr2:Segkk3|quiver*->Segkk19|quiver->Segkk8|quiver*`

In fact, the configuration string I used for *D. yakuba* Tai18E2 before Pilon polishing was:

`X:Segkk0_quiver*->Segkk61_quiver*->Segkk60_quiver->Segkk41_quiver*->Segkk49_quiver=>2L:Segkk46_quiver*->Segkk47_quiver=>2R:Segkk7_quiver*=>3L:Segkk1_quiver*->Segkk54_quiver=>3R:Segkk48_quiver->Segkk55_quiver=>4:Segkk3_quiver*`

### `prepContigsForGenBank.pl`

Usage:

`prepContigsForGenBank.pl `

This script prepares most of the files you'll need for tbl2asn to submit a genome to GenBank.  It takes in a FASTQ (make sure that the contig names comply with GenBank requirements), and outputs an AGP (with filename specified as an argument) and a renamed FASTQ (on STDOUT) in case you need to annotate/relabel any contigs in particular (e.g. mtDNA contig needs a special label to specify a different genetic code).

As seen above with `manualScaffold.pl`, you specify the scaffolding order for the AGP with a configuration string.  Another requirement for submitting to GenBank is that scaffolds cannot have the same name as a contig used to make them.  This is an AGP requirement (though I didn't find it in the format specification...).  Unscaffolded contigs will get placed in the AGP under a scaffold name that is simply the contig name with _scaf appended.

You can specify the gap length and gap type as well, although last I checked NCBI wants them to always be type U and length 100 if you don't actually have an estimate for the gap size.

### `configStringToAGP.pl`

Usage:

`configStringToAGP.pl -i [.fai file] [-c config string file] [-u] [-a output AGP] <config string>`

`-i` specifies the FASTA index (.fai) for the unscaffolded FASTA, which is necessary to detect contigs missing from the config string.

`-c` acts as an alternate to providing the config string as a positional argument.

`-u` ensures that contigs missing from the config string are included in the output AGP.

`-a` specifies the path for the output AGP. If omitted, this defaults to STDOUT.

As you can imagine, this script basically does one part of what `prepContigsForGenBank.pl` does.  In fact it was derived from the same code, I just didn't need a FASTQ, just an easy way to get an AGP for quick manual scaffolding.

Note the slight differences between an assembly scaffolded with an AGP produced by `configStringToAGP.pl` versus scaffolding with the config string: The AGP version 2.0 specification states that gaps of unknown size (`U` records) should have length 100, whereas `manualScaffold.pl` using a config string will by default make gaps of 500 Ns.

### `fastqToFastaQual.pl`

Usage:

`fastqToFastaQual.pl `

Just another script for GenBank submission, where tbl2asn does not accept FASTQs, so you need to generate a .fsa, and optionally a .qual.  This script takes the FASTQ generated by `prepContigsForGenBank.pl` and creates a .fsa and a .qual.  If you don't want to submit with quality scores, make sure to rename the .qual file (or delete it).

### `softmaskFromHardmask.cpp`

This program will softmask a genome, given an unmasked genome and a hardmasked genome. The underlying code is quite simple, but it turns out to be a lot easier to use this than to re-run RepeatMasker on softmasking mode if you accidentally ran hardmasking mode. Just make sure both the unmasked and hardmasked FASTAs are either fully unwrapped, or wrapped to exactly the same length, otherwise the program will error out. Also, it doesn't seem to handle process substitutions correctly (e.g. will error out about different line wrappings if you feed it process substitutions calling fasta_formatter).

Usage:

`softmaskFromHardmask [unmasked FASTA] [hardmasked FASTA]`

## Annotation-related scripts:

### `ScaffoldGFF3.pl`

Usage:

`ScaffoldGFF3.pl [-b <Output Broken GFF3>] -a <AGP> -g <GFF3> [-i]`

`-a` specifies an AGP file indicating the relationships between the scaffolds in the GFF3 and the chromosomes they compose.

`-g` specifies the GFF3 annotation in whichever space is input (the more-scaffolded space if `-i` is used, less-scaffolded space otherwise)

`-b` specifies the filename for the output GFF3 of features unable to be converted between coordinate spaces

`-i` is a flag meant for inverting the direction of the coordinate space transformation specified by the AGP, so if the AGP goes from scaffolds to chromosomes, using `-i` would take a chromosomal GFF3 as input, and produce a scaffold-space GFF3.

The coordinate-transformed GFF3 is output to STDOUT.

The `-d` flag can be used multiple times, and triggers output of debugging information onto STDERR.

Example Usage:

`ScaffoldGFF3.pl -dd -b Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_toArms_broken.gff3 -a Dsan_STOCAGO1482_chromosome_arms_nocontam_mtDNA.agp -g Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_contigs.gff3 > Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_arms.gff3`

### `ScaffoldGFFtoChromosomeGFF.pl` (deprecated in favor of `ScaffoldGFF3.pl`)

**NOTE: This script is deprecated in favor of using `ScaffoldGFF3.pl`, which has more functionality, and likely bug fixes.**

Usage:

`ScaffoldGFFtoChromosomeGFF.pl -a [AGP] -g [GFF3]`

`-a` specifies an AGP file indicating the relationships between the scaffolds in the GFF3 and the chromosomes they compose.

`-g` specifies the GFF3 annotation in scaffold-space that needs to be converted to chromosome-space.

Example Usage:

`ScaffoldGFFtoChromosomeGFF.pl -a Hmel2.agp -g Hmel2.gff3 > Hmel2_chromosomes.gff3`

This was just a quick script I wrote to go from an annotation (GFF) in scaffold coordinate space (in my case, it was the *Heliconius melpomene* v2 assembly) into chromosomal coordinates based on an AGP.  Not sure how many others will need to use this, but it came in handy for Hmel2.

### `extractRegions.pl`

Usage:

`extractRegions.pl -i [FASTA] -b [BED] [-p [PREFIX]]`

`-i` is optional and defaults to STDIN

`-p` is optional and defaults to no prefix for FASTA headers

Example Usage:

`extractRegions.pl -i Dyak_NY73PB_v2_w60.fasta -b Dyak_NY73PB_v2_4fold_sites_genomic.bed -p Dyak_NY73PB_v2 > Dyak_NY73PB_v2_4fold_sites.fasta`

This script extracts specific segments of a FASTA based on the intervals provided by a BED file.  If you extract only exon records from a GFF and convert that to BED, you could then extract all the exons as separate FASTA records.  Same thing with introns, assuming they have records in the GFF.  This is also useful for extracting flanking regions for primer design, etc.

### `extractSites.pl`

Usage:

`extractSites.pl -i [FASTA] -b [BED]`

`-i` is optional and defaults to STDIN

Example Usage:

`extractSites.pl -i Dyak_NY73PB_v2_w60.fasta -b Dyak_NY73PB_v2_4fold_sites_genomic.bed > Dyak_NY73PB_v2_4fold_sites.fasta`

The idea here is similar to `extractRegions.pl`, but this script is the predecessor, and was written out of necessity to create a concatenated FASTA record from many identified sites.  I used this script in combination with `identify4foldDegenerateSites.pl`, but that hasn't been uploaded yet.

### `overlappingFeatures.pl`

Usage:

`overlappingFeatures.pl `

Example Usage:

`overlappingFeatures.pl `

This script takes in a BED file of intervals you find interesting, and a GFF of features.  It then extracts lines of the GFF that overlap your intervals, and outputs a modified set of GFF lines.  I used this as a preliminary way to screen for genes within tracts of IBD common between multiple white monarchs.

### `fillMissingGenes.awk`

This awk script takes a GFF3 file produced by gtf2gff.pl (from Augustus, tested with Augustus version 3.3), and fills in any missing `gene` and `mRNA` records, which happens surprisingly often. This script has primarily been used in the process of converting the `augustus.hints.gtf` file produced by the BRAKER2 (v2.1.0) pipeline into a specification-compliant GFF3.

Example Usage:

`cat <(awk 'BEGIN{print "##gff-version 3";}{print "##sequence-region "$1" 1 "$2;}' Dsim_w501_Pilon_chromosome_arms_mtDNA_softmasked_w60.fasta.fai) <(fillMissingGenes.awk Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3) > Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted.gff3`

(Note that the first awk command in the example is simply adding in the appropriate sequence-region header lines based on the FASTA index file, as well as the GFF version header.)

The resultant GFF3 (with GFF version and sequence-region headers) should now be compatible with GenomeTools' sorting algorithm, like so:

`gt gff3 -sort -tidy -retainids -checkids Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted.gff3 2> Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted_sorting.stderr > Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted_sorted.gff3`

The command used to generate the `Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3` input file was:

`awk 'BEGIN{FS="\t";OFS="\t";}$2=="AUGUSTUS"{if ($3 == "transcript") {split($9, genetxarr, "."); $9="gene_id \""genetxarr[1]"\"; transcript_id \""$9"\";";}; print $0;}' braker/Dsim_w501_BRAKER2/augustus.hints.gtf | ~/bin/augustus-3.3/augustus/scripts/gtf2gff.pl --printExon --gff3 --out=Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3 2> Dsim_w501_BRAKER2_RNAseq_DmelProteins_gtf2gff.stderr`

### `amendTrinotateKEGG.pl`

This Perl script takes in a Trinotate XLS report file, and uses the REST API for KEGG to create a GOseq-compatible file of KO pathways for each transcript in the Trinotate report. If you already have the GOseq-compatible list of KOs for each transcript, you can use the `-g` argument to pass this in, and the script will output a gene-level version of the same file (combining the lists for alternative isoforms) without accessing the KEGG REST API.

Since this takes advantage of the KEGG REST API, please don't use it too often, as KEGG is kindly providing this resource to us, so we don't want to cause them trouble.

### `getKEGGdescriptions.pl`

This Perl script takes the output of `amendTrinotateKEGG.pl` and uses the KEGG REST API to produce a TSV of KO pathway IDs and their associated descriptions.

Since this takes advantage of the KEGG REST API, please don't use it too often, as KEGG is kindly providing this resource to us, so we don't want to cause them trouble.

### `constructCDSesFromGFF3.pl`

Usage:

`constructCDSesFromGFF3.pl -i [FASTA] -g [GFF3] [-e]`

`-i` specifies the FASTA from which to extract sequences. `-i` is optional, defaults to STDIN

`-g` specifies the GFF3, which must be in the same coordinate space as the FASTA, and must be locally sorted (i.e. exons and CDSes within a gene must be in ascending positional order).

`-e` is optional, outputs an extra 2 lines per FASTA record, which are the "Exon Range String" and the "CDS Range String".  These indicate the coordinates of the boundaries of each CDS record in genomic-space and CDS-space, respectively. They are always written with the initial exon first, and terminal exon last, so CDSes on the negative strand have Exon Range String elements listed in reverse order, and with larger number first in each range.

This script takes in a GFF3 file (e.g. genome annotation -- be very certain that it is GFF3 and not GFF2 or GTF), and a genome FASTA, and outputs all CDSes from the annotation as unwrapped FASTA.

I've used this to prepare transcriptomes (minus UTRs) for differential expression analysis, and it does not suffer from the off-by-one error you get from bedtools getfasta (plus you don't get all the extraneous individual exon and CDS FASTA records).

For instance, if you have the "GFF" output from the BRAKER1 pipeline, this is not actually GFF, but an oddly formatted GTF, so you can convert from GTF to GFF3 using the Augustus script `gtf2gff.pl` as follows:

```
gtf2gff.pl --printExon --printUTR --gff3 --out BRAKER_output_gtf2gff.gff3 < BRAKER_output.gff
constructCDSesFromGFF3.pl -i BRAKER_genome.fasta -g BRAKER_output_gtf2gff.gff3 > BRAKER_output_gtf2gff_transcripts.fasta
```

Oftentimes with BRAKER1- or BRAKER2-derived GFF3 files, you'll want to fix and fully sort them. My typical command set for getting a valid, sorted GFF3 ready for this script is:

```
awk 'BEGIN{FS="\t";OFS="\t";}$2=="AUGUSTUS"{if ($3 == "transcript") {split($9, genetxarr, "."); $9="gene_id \""genetxarr[1]"\"; transcript_id \""$9"\";";}; print $0;}' braker/[BRAKER run ID]/augustus.hints.gtf | ~/bin/augustus-3.3/augustus/scripts/gtf2gff.pl --printExon --gff3 --out=[annotation ID].gff3 2> [annotation ID]_gtf2gff.stderr
samtools faidx [softmasked genome FASTA]
cat <(awk 'BEGIN{print "##gff-version 3";}{print "##sequence-region "$1" 1 "$2;}' [softmasked genome FASTA].fai) <(fillMissingGenes.awk [annotation ID].gff3) > [annotation ID]_adjusted.gff3
gt gff3 -sort -tidy -retainids -checkids [annotation ID]_adjusted.gff3 2> [annotation ID]_adjusted_sorting.stderr > [annotation ID]_adjusted_sorted.gff3
```

### `extractIntronsFromGFF3.pl`

Usage:

`extractIntronsFromGFF3.pl -i [FASTA] -g [GFF3]`

`-i` specifies the FASTA from which to extract sequences. `-i` is optional, defaults to STDIN

`-g` specifies the GFF3, which must be in the same coordinate space as the FASTA, and must be locally sorted (i.e. exons and CDSes within a gene must be in ascending positional order).

This script takes in a GFF3 file (e.g. genome annotation -- be very certain that it is GFF3 and not GFF2 or GTF), and a genome FASTA, and outputs all CDSes from the annotation as unwrapped FASTA.

The FASTA headers are space-separated lists of information useful for downstream filtering for homology: Intron ID (GFF3-style), intron coordinates (FlyBase-style), intron length, left flanking exon length, right flanking exon length, left flanking exon ID, right flanking exon ID.

Beware that the extracted intronic sequence will contain the two splice junctions, so you'll probably want to trim those before aligning.

GFF3-style intron and exon IDs are of the form: [transcript ID].intron# or [transcript ID].exon#

FlyBase-style coordinates are of the form: `[scaffold ID]:[left bound]..[right bound]([strand])`

### `createGFFfromExonerate.pl`

Usage:

`createGFFfromExonerate.pl -v [Exonerate output with VULGAR strings] -f [FAI of genome] [-w wonky genes GFF3] [-p prefix]`

`-v` specifies the output from Exonerate when mapping proteins to a genome with `--model protein2genome`, specifically requiring `--showvulgar yes --showalignment yes`.

`-f` specifies the FASTA index (.fai) for the genome FASTA the proteins were mapped to. This file is used to generate a valid GFF3 header, but does not impact the contents of any of the GFF3 records.

`-w` specifies an output file for "wonky" genes which don't fit the typical model of a gene. TODO: Elaborate on what makes them wonky.

`-p` specifies the prefix used for gene, transcript, and protein IDs. For example, FlyBase labels their genes as FBgn######, their transcripts as FBtr######, and their proteins as FBpp######, so the prefix there is "FB". For the KAIKObase Bombyx mori annotation, the prefix is "Bm".

Example Usage:

`createGFFfromExonerate.pl -v DyakTai18E2_DmelProteinsExonerate.out -f Dyak_Tai18E2_renamed_w60.fasta.fai -w Dyak_Tai18E2_DmelProteins_wonky_genes.gff3 -p "FB" > Dyak_Tai18E2_DmelProteins_fromExonerate.gff3`

This script takes VULGAR alignments from an Exonerate run of proteins against a genome, and generates a [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) of coding sequences, including all gene and mRNA records, if the script is able to infer the appropriate gene and/or mRNA ID.  The gene ID is usually scraped from the human-readable alignment portion of Exonerate's output (in particular, from the contents of the "Query: " line), looking for a pair of IDs after a "parent=" tag definition, as is seen in FlyBase or KAIKObase headers.  The first ID is assumed to be the gene ID, and the second is assumed to be the transcript ID.  If no such parent= tag can be found, new gene and transcript IDs are synthesized by prefixing the protein ID (prior to the first space in its name) with a supplied prefix (default: my) and "gn" for gene IDs, and "tr" for transcript IDs.

Any protein alignments with particularly odd configurations won't be output to the main GFF3, and are instead optionally output to a file specified with the `-w` option.

### `extractCDSfromVulgar.pl`

Usage:

`extractCDSfromVulgar.pl -v [Exonerate output] -i [Genome FASTA] -p [prefix]`

`-v` specifies the output from Exonerate when mapping proteins to a genome with `--model protein2genome`, specifically requiring `--showvulgar yes --showalignment yes`.

`-i` specifies the genome FASTA, from which nucleotides are extracted to compose each CDS.

`-p` specifies a prefix to prepend to the header of each CDS.

Example Usage:

`extractCDSfromVulgar.pl -v DyakTai18E2_DmelProteinsExonerate.out -i Dyak_Tai18E2_renamed_w60.fasta -p Dyak_Tai18E2 > Dyak_Tai18E2_DmelProteins_CDSes.fasta`

This script takes VULGAR alignments from an Exonerate run of proteins against a genome, and extracts CDSes for each protein provided, skipping frameshifts (and implicitly deletions) so that the resulting CDS is still in-frame.

### `codingSitesByDegeneracy.pl`

Usage:

`codingSitesByDegeneracy.pl -f [Degree of Degeneracy] -i [CDS FASTA]`

`-f` specifies the degree of degeneracy of sites you want to identify. This degree is the number of nucleotides for which, when a nucleotide is substituted at this position, the amino acid encoded by the codon remains the same as the input codon. For example, at a 4-fold site (degree of degeneracy 4), any of the 4 nucleotides may be substituted at this position to obtain the same amino acid.

`-i` specifies the path to a FASTA of coding sequences (CDSes) whose sites you want to categorize.

Example Usage:

`codingSitesByDegeneracy.pl -f 4 -i Dyak_NY73PB_v2_BRAKER2_RNAseq_DmelProteins_CDSes.fasta | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > Dyak_NY73PB_v2_4fold_sites.bed`

This script takes a degree of degeneracy (specified by the `-f` flag, e.g. `-f 4` for 4-fold synonymous sites, `-f 2` for 2-fold sites, etc.) and a FASTA of coding sequences (CDSes, say generated by `constructCDSesFromGFF3.pl`) and outputs an unsorted and unmerged BED indicating the CDS-space positions of sites with the specified degree of degeneracy.

When coupled with `CDStoGenomicIntervals.pl`, `subsetVCFstats.pl`, and `calculateDxy`, it allows for calculation of an uncorrected per-site dS (or even dN for `-f 0` or `-f 1`).

### `CDStoGenomicIntervals.pl`

Usage:

`CDStoGenomicIntervals.pl -i [CDS-space BED] -g [GFF3]`

`-i` specifies a BED file of intervals in CDS-space which you want to convert to genome-space.

`-g` specifies the GFF3 file that provides the mapping from CDS-space to genome-space.

Example Usage:

`CDStoGenomicIntervals.pl -i Dyak_NY73PB_v2_4fold_sites.bed -g Dyak_NY73PB_v2_BRAKER2_RNAseq_DmelProteins_adjusted_sorted.gff | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > Dyak_NY73PB_v2_4fold_sites_genomic.bed`

This script takes a BED of CDS-space intervals and a GFF3 annotation (requiring CDS records that have a Parent= tag whose value would match the scaffold ID of the CDS-space BED), and outputs the equivalent BED in genome-space, accounting for intron-spanning intervals.

The output BED is **not** sorted or merged, so pass the output through `sort -k1,1 -k2,2n -k3,3n` and `bedtools merge -i -` before downstream usage.

## Genome-wide statistics programs/scripts:

**Note: All C++ programs will safely compile as long as your compiler supports C++11.  I usually compile with `g++ -O3 -g -Wall --std=c++11 -o [program prefix] [program prefix].cpp`.**

### `listPolyDivSites.cpp`

This program expects two FASTAs with identical lengths and identical line-wrap lengths.  It outputs a simple 3 or 4 column TSV, one line per base, with columns as follows:

1. Scaffold ID
2. Position in reference FASTA
3. 1 or 0 indicating whether or not the sample FASTA is polymorphic or divergent relative to the reference FASTA
4. 1 or 0 indicating whether or not either FASTA has an N at this position

A typical call might look like this:

`listPolyDivSites -p -n [reference FASTA] [sample FASTA] | nonOverlappingWindows -n -w [window size in bp] -o [sample prefix]_poly_w[window size in kb]kb.tsv`

### `oneSamplePolyDiv.sh`

This script is a wrapper for finding windowed heterozygosity and rate of fixed differences for one sample against the reference.  It basically just wraps the example line for `listPolyDivSites` in the case of "poly" being specified, and uses the `-d` flag and a different output TSV filename for the "div" case.

Example usages:
`oneSamplePolyDiv.sh poly [sample FASTA] [reference FASTA] [window size] > [sample prefix]_poly_w[window size in kb]kb.tsv`
`oneSamplePolyDiv.sh div [sample FASTA] [reference FASTA] [window size] > [sample prefix]_div_w[window size in kb]kb.tsv`

### `oneSamplePolyDivPlot.R`

This is a quick R script for making single-sample heterozygosity and rate of fixed difference plots for windows across the genome.  I think I hard-coded the specific set of scaffolds that it subsets out (i.e. the chromosome arms of the *Drosophila melanogaster* subgroup), so you may have to edit that for your own purposes.  And maybe having all scaffolds on the same plot won't suit your purposes.  But it shouldn't be hard to edit the ggplot call to fix that.

Example usage, given you ran oneSamplePolyDiv.sh as in the examples above:
`oneSamplePolyDivPlot.R ./ [sample prefix] [short ID for reference] [window size in bp] [maximum y coordinate to plot] arms`

This produces a PDF with heterozygosity and fixed difference rate across the Drosophila chromosome arms (excluding the 4). The PDF filename follows the template `[sample prefix]_polydiv_w[window size in kb]kb.pdf`.

Options for the last argument:

1. `arms`: equivalent to `X,2L,2R,3L,3R`
1. `allArms`: equivalent to `X,2L,2R,3L,3R,4`
1. `chroms`: equivalent to `X,2,3`
1. `allChroms`: equivalent to `X,2,3,4`
1. comma-separated list of scaffolds or chromosomes

### `oneSamplePolyDivPerScaf.sh`

Example Usage:

`oneSamplePolyDivPerScaf.sh poly CY01A_realigned_MPILEUP_final_pseudoref.fasta Dyak_NY73PB_v2_w60.fasta`
`oneSamplePolyDivPerScaf.sh div CY01A_realigned_MPILEUP_final_pseudoref.fasta Dyak_NY73PB_v2_w60.fasta`

Rather than finding windowed values for heterozygosity and rate of fixed differences, this wrapper averages the statistics across each scaffold.

### `oneSamplePolyDivOverall.sh`

Example Usage:

`oneSamplePolyDivOverall.sh poly CY01A_realigned_MPILEUP_final_pseudoref.fasta Dyak_NY73PB_v2_w60.fasta`
`oneSamplePolyDivOverall.sh div CY01A_realigned_MPILEUP_final_pseudoref.fasta Dyak_NY73PB_v2_w60.fasta`

Even broader, this calculates genome-wide heterozygosity and rate of fixed differences for a single sample against the reference.

### `nonOverlappingWindows.cpp`

**Version change:** As of version 1.2, nonOverlappingWindows automatically skips the first line if it is a header line (i.e. does not contain any numbers), and has an option to indicate which column of the file to use as a statistic column.

Usage:

`nonOverlappingWindows [-n] [-w window size in bp] [-i input TSV] [-o output TSV] [-s stat column]`

`-n` indicates whether or not the input TSV has a fourth column indicating whether or not to omit the current site (1=omit, 0=keep).

`-w` specifies the window size to use, in base pairs.

`-i` specifies the path to the input TSV, which consists of columns: Scaffold ID, Position, Statistic, and optionally an Omit column. This argument is optional, the default is STDIN.

`-o` specifies the output TSV, which consists of columns: Scaffold ID, Position, Averaged Statistic, and an optional column indicating the fraction of used sites (i.e. 1 - fraction of omitted sites).

`-u` triggers the output of the fourth column, i.e. fraction of usable sites

`-s` specifies the column to use for averaging the statistic. The default is column 3 (i.e. `3`), but can be any integer above 2 as long as it specifies a valid column.

Example Usage:

`listPolyDivSites -n -p Dyak_NY73PB_v2_w60.fasta CY01A.fasta | nonOverlappingWindows -n -w 100000 -o CY01A_poly_w100kb.tsv`

This program takes a TSV with 3 or 4 columns, and averages the values in the 3rd column over windows of specified length.  Using the `-n` flag leads to omission of positions within the window that have a 1 in the 4th column of the input.  These sites are omitted from both the numerator and denominator of the average, hence we don't bias the estimate by interpreting masked bases as anything other than missing data.

It was originally written to calculate windowed depth using the output of `samtools depth -aa`, but the input format is general enough that most if not all of my stats tools use it.

### `calculateDxy.cpp`

**Version change:** As of version 2.2, you do not need to list the FASTAs as positional arguments, as the paths to the FASTAs are read from the populations metadata file. This makes for a substantially shorter command line.

Among the many basic stats we might want to calculate, Dxy and Pi are pretty basic.  This program calculates both, given a TSV that maps FASTA filenames to population numbers, and a list of FASTA filenames as positional arguments. The output has a variable number of columns, dependent on the number of populations specified.  The first four columns will always be:

1. Scaffold ID
2. Position
3. Dxy between populations 1 and 2
4. Whether this site should be omitted (or, with the `-u`, is the fraction of usable sequences at that position)

The remaining columns are Dxy and Pi for the various combinations of populations.

Dxy is calculated on a per-site basis, using the generalized equation (pardon the lack of LaTeX interpretation here):

D\_{x, y} = \frac{n}{n-1} \sum\_{i,j \in {A,C,G,T}} 2\*p\_{x,i}\*p\_{y,j}

And Pi is calculated as:

Pi\_{x} = \frac{n}{n-1} \sum\_{i < j} 2\*p\_{x,i}\*p\_{x,j}

These are no longer exclusively for biallelic sites, but do reduce to the biallelic site formulae when i and j are constrained to 1 and 2.  Also, I haven't actually written out the proof yet, but I'm pretty sure that averaging the per-base values across a window is exactly equivalent to the Dxy and Pi you would calculate using Nei's formulae.

`n` in these formulae is the haploid sample size, and we assume that diploid sequences are provided (heterozygous sites identified by IUPAC degenerate bases K, M, R, S, W, and Y. Thus, without the `-i` flag, `n` is 2 times the number of input sequences.

The `-i` flag indicates that only 1 allele per site per input FASTA should be counted in the allele frequencies. The allele is chosen at random for each input sequence at each site (the `-r` argument sets the PRNG seed), and the `n` in the above formulae is equal to the number of input sequences. This is appropriate for calculations based on inbred or partially inbred individuals.

Note that you have to omit the first line (a header) in order to pass the output to `nonOverlappingWindows`.  Also, if you want to take the windowed average of something other than D12, you'll need to pipe the output through a quick awk script, e.g. to average windows of the 6th column:

Further note: Versions 2.2 and up of `nonOverlappingWindows` will automatically skip the above-mentioned header line, and have the `-s` argument for selecting which column has the statistic of interest.

`calculateDxy -p [population map TSV] | awk 'NR>1{print $1"\t"$2"\t"$6"\t"$4;}' | nonOverlappingWindows -n -w [window size in bp] -o [output TSV filename]`

or in newer versions:

`calculateDxy -p [population map TSV] | nonOverlappingWindows -s 6 -n -w [window size in bp] -o [output TSV filename]`

The usable fraction of sequences can be used as a filter for the windowed averaging, e.g. for a weighted average, using the usable fraction as the weight for a site (the `-a` flag for `nonOverlappingWindows`):

`calculateDxy -p [population map TSV] -u | nonOverlappingWindows -a -w [window size in bp] -o [output TSV filename]`

### `calculatePolymorphism.cpp`

This program calculates pi given a list of FASTA filenames as positional arguments. The output columns are:

1. Scaffold ID
2. Position
3. Pi within the samples
4. Whether this site should be omitted (or, with the `-u`, is the fraction of usable sequences at that position)

Pi is calculated as in calculateDxy.

These are no longer exclusively for biallelic sites, but do reduce to the biallelic site formulae when i and j are constrained to 1 and 2.  Also, I haven't actually written out the proof yet, but I'm pretty sure that averaging the per-base values across a window is exactly equivalent to the Dxy and Pi you would calculate using Nei's formulae.

`n` in these formulae is the haploid sample size, and we assume that diploid sequences are provided (heterozygous sites identified by IUPAC degenerate bases K, M, R, S, W, and Y. Thus, without the `-i` flag, `n` is 2 times the number of input sequences.

The `-i` flag indicates that only 1 allele per site per input FASTA should be counted in the allele frequencies. The allele is chosen at random for each input sequence at each site (the `-r` argument sets the PRNG seed), and the `n` in the above formulae is equal to the number of input sequences. This is appropriate for calculations based on inbred or partially inbred individuals.

The `-s` flag causes the program to output a 1 if the site is segregating (i.e. pi > 0.0), and a 0 if not. This is the simple logical extension of `listPolyDivSites -p` to multiple samples, and can be used to evaluate Watterson's estimate of theta (you can calculate the number of segregating sites, S, sometimes denoted k, quickly from the output).

Note that you have to omit the first line (a header) in order to pass the output to `nonOverlappingWindows`, e.g. using `tail -n+2` or awk, as follows:

Further note: Versions 2.2 and up of `nonOverlappingWindows` will automatically skip the above-mentioned header line, and have the `-s` argument for selecting which column has the statistic of interest.

`calculatePolymorphism [FASTA 1] [FASTA 2] [FASTA 3] [...] | awk 'NR>1' | nonOverlappingWindows -n -w [window size in bp] -o [output TSV filename]`

or with newer versions of `nonOverlappingWindows`:

`calculatePolymorphism [FASTA 1] [FASTA 2] [FASTA 3] [...] | nonOverlappingWindows -n -w [window size in bp] -o [output TSV filename]`

In the degenerate case of inputting a single FASTA, this program behaves like `listPolyDivSites -p -n`, outputting 1 for heterozygous sites, 0 for all others.

### `subsetVCFstats.pl`

Usage:

`subsetVCFstats.pl -i [input TSV] -b [BED]`

`-i` specifies the path to the input TSV file (where the first two columns are Scaffold ID, and Position). Lines beginning with `#` will be skipped. If omitted, this argument defaults to STDIN.

`-b` specifies the BED file containing intervals to keep/subset from the input TSV.

Example Usage:

`subsetVCFstats.pl -i Dyak_Dsan_Dxy.tsv -b Dyak_4fold_sites_genomic.bed > Dyak_Dsan_Dxy_4fold_sites.tsv`

This script takes an input TSV file (where the first two columns are scaffold ID and position), and a BED file, and subsets out lines of the TSV that correspond to intervals in the BED file.  It is generally useful for subsetting, whether subsetting lines from a VCF, or Dxy or pi values from the output of `calculateDxy`, or per-base depths from the output of `samtools depth`, etc.

### `decompressStats.pl`

This script adds in missing records into a subsetted statistics TSV so that the output can be used by `nonOverlappingWindows`. Any sites added in have the Omit column set to 1 so that they are not counted in the average for a window.

Usage:

`decompressStats.pl [-i stats TSV] [-f genome .fai] [-s statistic column]`

`-i` specifies the input per-site statistics TSV, with columns being Scaffold ID, Position, Statistic, and Omit. If omitted, the default is STDIN.

`-f` specifies the path to the FASTA index (.fai generated by `samtools faidx`) for the original genome.

`-s` specifies the column to use as a statistic (default is 3, but can otherwise be any integer greater than 4).

Example Usage:

`subsetVCFstats.pl -i Dyak_genomewide_Dxy.tsv -b Dyak_4fold_sites.bed | decompressStats.pl -s 3 -f Dyak.fai > Dyak_4fold_Dxy.tsv`

### `divergenceFromMAF.pl`

This script calculates per-site haploid divergence between two references aligned by whole-genome alignment, based on an alignment of 1:1 regions in MAF format. For example, a MAF generated by the LAST protocol used for Human-Chimp alignments, filtering for 1:1 alignments, would be appropriate as input.

Usage:

`divergenceFromMAF.pl [-i input MAF] [-s] <Species 1 Prefix> <Species 2 Prefix>`

`-i` specifies the path to the input MAF alignment file. Default is STDIN.

`-s` indicates whether the output should be in sorted order.

`Species 1 Prefix` is the prefix prepended to scaffold IDs of Species 1. The expectation here is that Scaffold IDs prior to alignment (or perhaps during identification of 1:1 alignments) follow the pattern `[Prefix].[ID]`, e.g. `DyakTai18E2.X` or `DteiGT53w.3R`. Thus, in these examples, you would specify `DyakTai18E2` or `DteiGT53w` as the species prefixes, respectively.

`Species 2 Prefix` is defined analogously to `Species 1 Prefix`, it just defines the pair in the case of `n > 2` whole-genome alignments.

Which species you choose to be `Species 1` determines the coordinate space of the output. That is, if you align `DyakTai18E2` and `DyakNY73PB`, and you choose `DyakTai18E2` as `Species 1 Prefix`, the output will be in `DyakTai18E2` coordinate space.

The contiguity of the output depends on the contiguity of the alignment blocks, so you will likely have to use `decompressStats.pl` before passing the output to `nonOverlappingWindows` or similar windowing scripts.

Example Usage:

`divergenceFromMAF.pl `

**Documentation of divergenceFromMAF.pl to be completed**

### `sitePatterns.cpp`

This program outputs counts of the genotype configurations found in a given set of FASTAs in the same coordinate space. Its output can be further processed for purposes like determining the empirical frequency of specific allelic configurations in an ASE experiment.

Usage:

`sitePatterns [list of pseudoreference FASTAs]`

## MSA-related scripts:

These scripts are related to pre-processing (and post-processing) multiple sequence alignments of coding sequences. They have been used as part of a pipeline to generate MSAs of about 9,400 single-copy orthologs across Dmel, Dsan, Dtei, and Dyak, and to integrate population resequencing data into these alignments.

### `addAlignedGaps.pl`

This Perl script takes in an aligned set of CDSes from reference genomes (aligned by, say, PRANK, and provided in aligned multi-FASTA format via the `-a` argument), and an unaligned multi-FASTA of population sample CDSes (provided with the `-i` argument, or on STDIN), matches up the species prefix in their FASTA headers (e.g. `>Dsan_` in the ref alignment would match to `>Dsan_` in the population samples), and inserts gaps into the population sample sequences to align them. We perform no optimization of the MSA, simply assume that the alignment of the references is correct, and add the population samples to the output alignment (file specified by flag `-o`, or STDOUT by default). This is FAR FAR faster than inserting all the unaligned reference and population sample sequences into an MSA program, without losing much information, especially for relatively low distance alignments.

Example Usage:

`addAlignedGaps.pl -a trimmed_aligned_refs/OG0001492_prank.best.fas -i trimmed_unaligned_samples/OG0001492_trimmed.fasta -o trimmed_aligned_all/OG0001492_untrimmed.fasta`

Note: I usually make sure the input FASTAs are unwrapped, but the code should work for arbitrary line-wrapping.  I just haven't tested it thoroughly.

### `fakeHaplotype.pl`

This Perl script takes in a diploid pseudoreference FASTA (e.g. produced by the [PseudoreferencePipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline/)) and randomly selects alleles at heterozygous sites to produce a haploid version of the input. The pseudo-random number generator (PRNG) seed can be set with the `-s` option (default: 42, because starting with the answer is better), and the `-a` flag allows you to choose the opposite allele, so you can run this script twice to produce the complementary two possible haploids from the input diploid. If you are splitting into both haplotypes, it is recommended to use the `-b` flag so that you can distinguish the fake haplotypes were you to concatenate the two outputs. The `-b` flag appends a `_0` or `_1` to the header of each sequence, `_0` if the `-a` flag is absent, `_1` if the `-a` flag is present. So this allows for safe splitting of haplotypes for non-haplotype-aware statistics from an outbred individual (or a synthetic diploid), and random selection of a haplotype from inbred individuals (make sure to only use one of the haplotypes, since the individual's ploidy is effectively less than 2).

Example Usage:

`cat <(fakeHaplotype.pl -i raw_CDSes_samples/OG0001492_unwrapped.fasta -b) <(fakeHaplotype.pl -i raw_CDSes_samples/OG0001492_unwrapped.fasta -a -b) | trimCDSes.awk > trimmed_unaligned_samples/OG0001492_trimmed.fasta`

This script is agnostic of line-wrapping, and will output in the same wrapping as the input.

### `skipRogers.awk`

This awk script is very very custom, and all it does is omit one of the two fake haplotypes from the Rogers et al. (2014) MBE Dyak data from the input alignment, and outputs all other sequences.

### `trimCDSes.awk`

This awk script takes in an unaligned set of CDSes (for example, extracted from reference assemblies, and organized into a single FASTA file per single-copy orthogroup, as identified by OrthoFinder analysis), trims them down to length that is an integral multiple of 3, and then removes a single terminal (i.e. 3'-most) stop codon, if detected, from each sequence. At the moment, only the stop codons in the NCBI transl_table=1 genetic code are detected.

Example Usage:

`trimCDSes.awk raw_CDSes_refs/OG0001492_unwrapped.fasta > trimmed_unaligned_refs/OG0001492_trimmed.fasta`

The script is extremely short and quick, enabling extreme parallelization with GNU Parallel.

### `trimAlignment.awk`

This awk script takes in an aligned FASTA, and truncates all sequences down to the same length as the shortest sequence in the set. For the output of a typical MSA program (e.g. Clustal Omega, MUSCLE, MAFFT, PRANK), this should amount to not doing anything to the input, but this script is useful when you use an alignment of reference sequences as a guide for inserting population samples into the alignment, and the population samples happen to have differing lengths. In particular, I wrote it for use with `addAlignedGaps.pl`, where a masked pseudoreference trimmed by `trimCDSes.awk` may not truncate the terminal stop codon (because `NNN` does not match `/T(AA|AG|GA)$/`), so some of the inserted pseudoreferences have 3 extra bases. By using this script, we rectify this problem in a very very simplistic manner.

Example Usage:

`trimAlignment.awk trimmed_aligned_all/OG0001492_untrimmed.fasta > trimmed_aligned_all/OG0001492_trimmed.fasta`

## Illumina data parsing scripts:

### `ReadHistogram.sh`

This quick bash script just simply wraps a one-liner that creates a histogram TSV of the index reads, sorted in descending order by count.  This histogram file is useful as a QC check before parsing, as it usually takes substantially less time to run than parsing.

It's essentially just two columns:

1. Index sequence
2. Count of index reads matching the sequence in column 1

Example usage:

`ReadHistogram.sh MySequencingLane_I1.fastq.gz > MySequencingLane_i7_index_histogram.tsv`

### `labelIndexReadHistogram.pl`

This script is intended to be used in tandem with `ReadHistogram.sh` (and may be piped to), in order to do a quick sanity check on your barcode file.  If the output from this script does not show most or all of your main barcodes with 0 mismatches at the top of the file, you either misspecified your barcodes file (maybe forgot to revcomp the index sequences?) or something went seriously wrong with the sequencing run.  This also gives a first-pass idea of how much error there was in the index read, and how many reads you should expect after parsing.

The output is a modified version of what comes from `ReadHistogram.sh` in that a variable number of columns is added, one per matching barcode from your barcodes file (with up to 2 mismatches).

Example usage:

`ReadHistogram.sh MySequencingLane_I1.fastq.gz | labelIndexReadHistogram.pl -b MySequencingLane_i7_indices.tsv > MySequencingLane_labeled_i7_histogram.tsv`

You may want to save the read histogram to a file, and feed it to `labelIndexReadHistogram.pl` with the `-i` option, in case you need to diagnose barcode file problems without waiting a long time for each run of `ReadHistogram.sh`.

**Note that your barcode file must be a proper TSV (i.e. columns separated by single `\t` characters, NOT spaces.  Many text editors have the annoying behaviour of inputting a certain number of spaces instead of a true tab character.  Using spaces will make the barcode parser output gibberish/fail, and will not produce any labels from this script.**

An easy way to check for tabs in your barcode file is to run `hexdump -C < MySequencingLane_i7_indices.tsv` and examine the output for `09` characters.  The ASCII code for the tab character is 09 in hexadecimal (see [ASCII Table](http://asciitable.com/)

### `divideConquerParser.sh`

This bash script was a simple (read: rushed) attempt to parallelize the slow process of parsing multiple plates out of an Illumina HiSeq 4000 lane.  At one point I was getting parsing jobs timing out after 24 hours, which seemed ridiculous given the embarassingly parallel nature of parsing.

All that happens behind the scenes here is that we split the read files into `n` parts, run `n` parallel instances of the barcode parser `barcode_splitter.py`, and once all are done, we re-merge the split files appropriately, and delete intermediate files.

Example usage:

Use 8 cores for a single-indexed paired-end (so R1, R2, and I1 FASTQ files) dataset:

`divideConquerParser.sh 3 "MyLane_R1.fastq.gz MyLane_R2.fastq.gz MyLane_I1.fastq.gz" 8 MyLane_i7_indices.tsv 3`

Use 8 cores for single-indexed paired-end dataset, but changing the listing order of the files so the index read file is first:

`divideConquerParser.sh 3 "MyLane_I1.fastq.gz MyLane_R1.fastq.gz MyLane_R2.fastq.gz" 8 MyLane_i7_indices.tsv 1`

Use 8 cores for single-indexed single-end dataset, where index read file is last:

`divideConquerParser.sh 2 "MyLane_R1.fastq.gz MyLane_I1.fastq.gz" 8 MyLane_i7_indices.tsv 2`

Use 8 cores to perform the i5 parse of a dual-indexed paired-end dataset (where i5 is the *_I1.fastq.gz file, and i7 is the *_I2.fastq.gz file):

`divideConquerParser.sh 4 "MyLane_R1.fastq.gz MyLane_R2.fastq.gz MyLane_I1.fastq.gz MyLane_I2.fastq.gz" 8 MyLane_i5_indices.tsv 3`

**Note that the 5th (last) argument is the 1-based position of the index read file you want to parse on.  This number must be less than or equal to the total number of read files you would like to parse.
