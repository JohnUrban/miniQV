![logo](/img/logo.png)

# miniQV
Assembly QV metrics based on fast alignments of a representative sample of reads.



# Detailed Usage

	MiniQV 0.0.2

        This produces a variety of alignment-based statistics that can be used to compare a set of assemblies.
        - Simply use the same read set and MiniQV parameters on each assembly.
        - Ultimately, the QV aspect of miniQV was designed for high accuracy short reads.
		- High accuracy contigs or high acc long reads should be fine).
                - Can technically be used on/for any alignments though: short reads, long reads, contigs, etc.
                - The same types of outputs will be produced, but interpretations will have to change accordingly.
                - Although error-prone reads can produce results that can be compared between assemblies, they may not be meaningful in an absolute sense.
                        - QV values will appear quite low even when correcting for the error-rate.
                - Designed with single-end reads (or single long reads) in mind.
                - Paired reads will work fine, but they return joint percent ID (and error rate) estimates based on both reads.
                        - The number of split alignments will be 2 for pairs that both map and don't split.
	- Nevertheless, the mis-assembly aspects of this pipeline were designed with long reads in mind (although short reads can work too).
		- Prevalence of mis-assemblies assessed by comparing the prevalence of split alignments among assemblies.
		- These comparisons are not totally released with version 0.0.1 though.



	Usage: 
		miniQV {-a asm.fasta | -b sam/bam } [-f] [-r READS.fastq(.gz)] [other options]

	Example Usage:

        If you have the target FASTA (e.g. an assembly) and a reads file in FASTQ format:
                miniQV -a /path/to/assembly.fasta -R /path/to/Reads.fastq.gz -s 1000000 -n 1000000 


        If you have a FOFN of multiple assemblies and a reads file (see -f notes below for more info on creating FOFN):
                miniQV -f -a /path/to/assembly.fofn -R /path/to/Reads.fastq.gz -s 1000000 -n 1000000 


        If you already have the alignments (any type*) in BAM format (from ANY aligner**):
                miniQV -b /path/to/alignments.bam
                *short reads, paired, single-end, long reads, contigs, etc.
                **Minimap2, BWA, bowtie2, Magic BLAST, STAR, HiSat2, etc


        If you already have an FOFN describing multiple BAMs/SAMs :
                miniQV -f -b /path/to/bams.fofn


        If you have an assembly FASTA and reads, and know the heterozygosity rate is 0.0204 (e.g. from GenomeScope2 report):
                miniQV -a /path/to/assembly.fasta -R /path/to/Reads.fastq.gz -s 1000000 -n 1000000 -H 0.0204


        If you have an assembly FASTA and reads, and know the heterozygosity rate is 0.0204, and want to save the SAM records:
                miniQV -a /path/to/assembly.fasta -R /path/to/Reads.fastq.gz -s 1000000 -n 1000000 -H 0.0204 -W intermediate.sam


        If you have an assembly FASTA and reads, and know the heterozygosity rate is 0.0204 (e.g. from GenomeScope2 report), and you are on SLURM:
                miniQV -a /path/to/assembly.fasta -R /path/to/Reads.fastq.gz -s 1000000 -n 1000000 -H 0.0204 -S -M 10G -T 10:00
	

        If you already have the alignments in BAM format, but they are position-sorted:
                miniQV -b /path/to/alignments.bam -N



	Options:
	
	-a	ASM	: Either -a or -b required. Provide /path/to/assembly.fasta OR /path/to/assemblies.fofn with -f to process multiple assemblies. See notes in -f.
			  See notes on file paths below.

	-b	BAM|SAM	: Either -a or -b required. 
			  Provide /path/to/readsPreAlignedToAssembly.bam[sam] 
				OR /path/to/bams.fofn with -f to process multiple BAM files. 
				See notes in -f.
			  IMPORTANT: For this analysis, reads need to be grouped or sorted by name, not position.
			     Most read mappers output all SAM/BAM lines for a given read together. 
				i.e. they are grouped by name (not sorted by name) and the output can be used directly with MiniQV.
			     If the BAMs are currently sorted by position, then use 'samtools sort -n pos-sorted.bam > name-sorted.bam'
				You can also use -N in miniQV to force name-sorting.
	    	             MiniQV allows taking just a sample of the lines from a BAM.
			     At the moment, sampling is simply taking a block of lines after optionally skipping a block of lines.
			     When the BAM is directly from the mapper, this block is essentially a random sample,
				 or at least it tends to give representative results compared to using all reads.
			     Whether or not that is true when the lines are sorted (by name) remains to be tested (TODO).
			     If using all the entire BAM file, or the majority of it, none of this matters: it will be representative.

        -f      FOFNMODE: Optional. No argument.
			  This flag lets miniQV know that the file provided to -a or -b is a special 2-col or 3-col file-of-filenames (FOFN).
			  Whether it is the 2-col or 3-col version is auto-detected by miniQV.
			  The 2-col version has columns:
				1 = prefix or nickname;
					used to make a sub-directory and as part of filenames for this assembly or bam.
				2 = absolute path to the asm.fasta or reads.bam file.
			  The 3-col version has columns:
			  	1 = group name ; all assemblies/bams with the same group name will be analyzed under that groupname subdirectory.
				    A use case is if you want to keep assemblies grouped together according to assembler used.
				    This directory structure is used with helper scripts to give the best assembly for each group (assembler), for example.
				2 = prefix or nickname; 
					used to	make a sub-directory within the group subdir and as part of filenames for this assembly or bam.
				3 = absolute path to the asm.fasta or reads.bam	file.

	-r	READS	: Required with -a. /Path/to/reads.fastq(.gz). 
			  Paired short reads can be analyzed separately, or concatenated into a single file.
			  In the latter case, the optimal number of split reads reported for a read name is 2 instead of 1. Here, 1 would indicate an orphaned read that had no splits.

	-d	DATATYPE: Possible values are 'sr', 'pacbio', and 'nanopore'. Default = sr (short read, presumably high accuracy).

	-m	MAPPER	: Possible values are 'minimap2', 'bowtie2', 'bwa', and 'magic'. 
			  Default = minimap2.
			  MiniQV was designed to be used with Minimap2 for all data types as part of a strategy to very quickly return QV estimates.
			  This is in part b/c Minimap2 can speedily generate indexes on the fly, and is an ultra fast mapper.
			  When minimap2 is combined with sampling just 1 million reads (for example), miniQV can generate evaluation results in minutes.
			  Using bowtie2/bwa will require relatively long assembly indexing steps in -a mode. 
			  However, if BAMs are already generated, then BAM mode (-b) skips mapping altogether, so it doesn't matter which aligner was used. 
			  Nevertheless, testing shows highly similar results between mappers. There is not a compelling reason to not use Minimap2, even on short reads

	-x	BTMODE	: Possible values are 'e2e' and 'local'. If using Bowtie2 as a mapper, then this tells it to use end-to-end or local mode.
			  Default = e2e (end to end mode).			  

	-t	THREADS	: Default = 16.

	-H	HETRATE	: Expected heterozygosity rate.
			  QV = -log10( error rate )
			  In aligned reads, the error rate is a sum of the errors from the assembly, the reads, and the natural variation (heterozygosity).
			  MiniQV attempts to show QVs with the total error rate as well as error rates after subtracting out error due to reads and/or natural variation.
			  It will show some results over a sweep of possible read error and heterozygosity rates in one output file.
			  However, for the cigar-QV table, it will only use the heterozygosity rate provided here.
			  Default = 0.001.
			  To get an expected heterozygosity rate, one method is to use GenomeScope2.

			  Coming soon -- some tools here to run that analysis automatically if needed.

	-s	SKIP	: Skip this number of reads in the beginning of the READS or BAM file.
                          For FASTQ files from Illumina, skipping at least the first million reads seems to be a best practice.
			    The error rate in the first million tends to be much higher than the average.
                          Only consider skipping reads if taken a sample (e.g. 1 million reads) with -n.
			  If using all of the reads (e.g. >20 million), no worries about skipping the first million.
			  Default = 0.

	-n	NREADS	: Only analyze the first N reads (after skipping reads if applicable).
			  With just 1 million reads, results are consistently highly representative of what the results would be when using all of the reads.
			    If just generating scores to compare assemblies, then the scores only matter relatively to each other anyway.
			    Thus, using the same block of 1 million reads on all assemblies allows meaningful direct comparisons very quickly.
			    Nonetheless, as said above, the results tend to be similar to using all reads when:
				 skipping the first million in a file, 
				 and then using the second million (-s 1000000 -n 1000000).
			  Default = 'all'.

	-S	SLURM	: Use to launch SLURM batch scripts (one per assembly in FOFN) instead of processing serially in current env. 
			  Default = false.

	-T	TIME	: Time limit for SLURM. Default = 12:00:00 .

	-M	MEM	: Mem limit for SLURM. Default = 80G .

	-W	SAVESAM : In -a mode, save the intermediate SAM alignments from the mapping step. 
			  Default = false (SAM records are just streamed into the next program).

	-N	SORT	: Force name sorting. 
			  This is only necessary if :
			  - using a mapper that doesn't output alignments grouped by read name, or 
			  - working with BAMs that are position sorted.

	-o	OUTDIR	: Name of main top-level output directory. 
			  Default = miniqv_output.

	-v	VERBOSE	: Default = false.





	Outputs:
	- See the *.err.txt and *out.txt files.
		- PRE.out.txt is a file that describes the alignment of every read. 
			- This is the output of samPctIdentity.py.
			- It has the number of splits, matches, mismatches, insertions, deletions, unaligned based.
				- these are used to create the *.cigar*.txt files.
			- It has 5 ways of computing percent identity. 
			- See samPctIdentity.py for more information.
		- miniQV_PRE.out.txt summarizes a few QV analyses over the entire assembly (see *cigar-QV.txt for even more depth, and contig info).
			- It shows estimates of error rates in the reads given their base qualities.
			- It summarizes how many reads aligned and various error sums.
			- It gives a completeness estimate here based on % of reads that aligned.
			- It gives various QV estimates, sweeping over an array of possible read error rates and heterozygosity rates:
				QV Estimates :
				- *Estimated QV from error rate in alignment blocks.
				- *Adjusted QV Estimates from error rate in alignment blocks plus unaligned bases (from aligned reads).
				- Adjusted QV Estimates from error rate in alignment blocks plus unaligned bases (from aligned reads) plus bases from unmapped reads.
				
				* found in cigar-QV.txt file too, using provided het rate and estimated read error rate

				Other estimates in this file that may be interesting, but are not the best QV estimates.
				- Estimated QV from proportion of imperfectly aligned reads as estimate of error rate.
				- Adjusted QV Estimates using proportion of imperfectly aligned reads as estimate of error rate given unaligned reads as well.
				- Percent Identity Method from Column 9 - Estimated QV from median error rate in aligned reads.
				- Percent Identity Method from Column 9 - Adjusted QV Estimates using median error rate in aligned given unaligned as well.
				- Percent Identity Method from Column 11 - Estimated QV from median error rate in aligned reads.
				- Percent Identity Method from Column 11 - Adjusted QV Estimates using median error rate in aligned given unaligned as well.
	- See the *cigar-QV.txt file.
		- This gives a table of QVs computed from various error rates (e.g. total error;   total error - read error;   total error - read error - hetrate)
		- It gives these QVs for entire assembly and for each contig.
	- See the stats files that give statistics on the percent identities of read alignments (min, max, mean, median, quantiles, etc).
		- these stats are the output of stats.py
		- miniQV automatically computes stats on columns 9 and 11 of the PRE.out.txt file.
	- See PRE.read-error.txt 
		- this information is also found in the miniQV.PRE.out.txt summary.	


	Use cases:

	- Created several assemblies ; want some metrics to compare them.
		- Just need a relative score.
		- The unadjusted QVs from the total assembly are fine for comparison.
			- No need to adjust for read error and het rate since those are constant.
			- Use same reads to analyze all assemblies (so same read error).
			- Assumes assemblies all have same underlying het rate.
				- e.g. all assemblies were made from same data set (different assemblers, parameters, polishers, etc).
				- e.g. all assemblies were made from same individual or population.
		- The completeness estimate based on % reads aligned can be used as well.


	- Have a final assembly, want an estimated QV value for the assembly (Assembly QV) and other absolute metrics.
		- Here you are looking for an absolute score.
		- The most conservative QV estimate is assuming 0 read error and 0 variation.
			- That is the unadjusted QVs.
		- Theoretically, the closest estimate to the true Assembly QV adjusts for read error and natural variation (heterozygosity).
			- miniQV attempts to do this by subtracting out :
				- the expected read error (given base call quality scores), and 
				- expected heterozygosity rate (e.g. from GenomeScope2 or other previous estimate).
			- miniQV outputs a maximum QV of 60 for now (1 error per million bases).
				- So if the QV after adjustment >60, it returns 60.
				- TODO: allow turning this feature off, or setting a different max.
		- The completeness estimate can be used as well.

	Citation:
	- This script was developed to evaluate the quality of various genome assemblies for the fungus gnat, Bradysia coprophila.
	- It is free to use, further develop, criticize, etc.
	- Feel free to cite this github repo if you find it helpful in your own research.



	Known issues:
	- File paths
			  miniQV should be installed inside a directory that has no spaces inside its absolute path.
			  In addition to miniQV, input file paths (-a, -b, -r) cannot have spaces. 
			  E.g. invalid paths = 
					/Path to my/file.fasta
					/Path\ to\ my/file.fasta
			       valid paths =
					/Path/to/my/file.fasta 
			  Spaces in paths are rare. They're more likely to be found on MacOS, and particular in Google Drive dirs.
			  Creating symlinks (ln -s) to the dir inside a path with no spaces can sometimes get around this.

			  The absolute path to miniQV should be used to execute miniQV (this is default if it is in your PATH). 
				Relative paths to miniQV can cause issues.
			  Files provided to -a, -b, and -r can be absolute paths or relative paths from the directory miniQV is executed inside.
			  	Still absolute paths are better than relative paths.
			  File paths inside the FOFN files should all be absolute paths (no relative paths).

	- Not enough memory (especially on SLURM).
		- If there is not enough memory available, the mapping stage of the pipeline can get killed.
		- This will lead to weird errors and outputs.
		- See ./known_issues/ for more info.



# Dependencies
- Minimap2
- Seqtk
- Samtools
- Python 2/3

# Other Dependencies (if going beyond defaults)
- Bowtie2
- BWA
- Magic BLAST




