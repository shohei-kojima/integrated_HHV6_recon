# This is a README file for iciHHV6_reconstruction_pipeline

# Flowchart
<img src='https://github.com/shohei-kojima/iciHHV6_reconstruction/blob/master/lib/image_for_readme.png' width='320px'>


# 0. prerequisites
### required software for running with default settings
- Linux (recommended: Ubuntu 18.04)

- Python 3.7 or later
- pysam 0.15.2 or later
- matplotlib 3.1.1 or later

- hisat2 2.1.0 or later 
- samtools 1.9 or later
- bcftools 1.9 or later
- bamCoverage 3.3.1 or later (deeptools 3.3.1 or later)
- gatk 4.1.7.0 or later
- picard 2.21.9 or later
- java (version matches to your picard)

### optional prerequisites
- bwa 0.7.17 or later (if you specify '-bwa' option)
- SPAdes genome assembler v3.13.1 or later (if you specify '-denovo' option)

### required Python built-in modules
os,sys,datetime,multiprocessing,logging,traceback,argparse,glob,pylab,subprocess,gzip



# 1. download this tool from GitHub (currently not public)
git clone https://github.com/shohei-kojima/iciHHV6_reconstruction



# 2. quick usage guide for the impatient

### when you use your BAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
In this case, you need to specify your BAM file with '-b' option. You also need to specify '-alignmentin' option as well.


### when you use your CRAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-c test.cram \
-fa /path/to/reference/genome/of/cram/GRCh38.fa \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
In this case, you need to specify your CRAM file with '-c' option. You also need to specify '-alignmentin' option and '-fa' option as well.

### when you use all discordant reads in BAM/CRAM file for mapping to the virus genome (alignmentin option + all_discordant option)
```
python main.py \
-alignmentin \
-all_discordant \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
By default, only unmapped reads for mapping to viruses when BAM or CRAM file is specified as an input. If you want to use all discordant reads (e.g. reads without sam flag '2') for mapping to viruses, you can specify '-all_discordant' option. The discordant reads includes read-pairs with distant mapped positions and ones with low MAPQ; therefore, the number of discordant reads are far higher than unmapped reads. This option is only available when '-alignmentin' option was specified. Although this option is less accurate than default setting (= use only unmapped reads), the coverage of reconstructed sequence is longer than default. Please use this at your own risk.

### when you use your fastq files as inputs (fastqin option)
```
python main.py \
-fastqin \
-fq1 test_1.fastq \
-fq2 test_2.fastq \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
You can specify your fastq files as inputs by '-fastqin' option. This tool does NOT support single-end reads. Please specify paired-end reads as separated files ('-fq1' and '-fq2' options).

### when you perform de novo assembly using metaspades (denovo option)
```
python main.py \
-alignmentin \
-denovo \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
'-denovo' option performs de novo assembly of HHV-6 sequence, in addition to variant call and consensus sequence generation. 



# 3. output files
### 'virus_detection_summary.txt'
This is one of the main result files for most users. This contains read coverage information of each virus genome.

- 1st column: RefSeq ID (fasta header without attribution)
- 2nd column: Whether a virus-derived reads are abundantly detected or not
- 3rd column: Coverage information
    - genome_length: Length of the virus genome
    - mapped_length: Length of the virus genome with one or more reads
    - perc_genome_mapped: Percent of virus genome with one or more reads
    - average_depth: Average of mapped read depth (average of whole genome)
    - average_depth_of_mapped_region: Average of mapped read depth of mapped regions (average of only mapped regions)
    - ratio_ave_virus_depth_to_autosome_depth: Average of mapped read depth of mapped regions devided by autosome depth provided with '-depth' option. Available only when '-depth' option was specified. Otherwise, 'NA'.
- 4th column: attribution of fasta header

### 'hhv6a_reconstructed.fa', 'hhv6b_reconstructed.fa'
This is one of the main result files for most users. This file contains HHV-6 sequence reconstructed with called variants. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B. Genomic regions where do not have any reads (= 0 reads mapped) are masked by a character 'N.'

### 'high_coverage_viruses.pdf'
This is one of the main result files for most users. If there is one or more viruses in your sample, this tool outputs read coverage of those viruses.

### 'hhv6a.vcf.gz', 'hhv6b.vcf.gz'
This is one of the main result files for most users. This file contains variant calls. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B.

### 'mapped_to_virus_dedup.bam'
BAM file containing alignment with virus genomes.

### 'mapped_to_virus.bedgraph'
Read depth of virus genomes.

### 'mark_duplicate_metrix.txt'
Summary of picard MarkDuplicates running with 'mapped_to_virus_dedup.bam.'

### 'hhv6a_metaspades_assembly', 'hhv6b_metaspades_assembly'
This file contains HHV-6 sequence reconstructed by metaspades. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B. You need to specify '-denovo' option to obtain this result.

### 'for_debug.log'
Log file. Stores all information, including input arguments and parameter setting.



# 4. options
### '-alignmentin'
Specify when you use BAM or CRAM file as an input. You also need to specify your input file with '-b' or '-c' option. When you use CRAM file as an input, you also need to specify reference genome with '-fa' option.

### '-b [BAM file]'
Use when specifing BAM file. Available only when '-alignmentin' is also specified.

### '-c [CRAM file] -fa [reference fasta file]'
Use when specifing CRAM file. Available only when '-alignmentin' is also specified.

### '-all_discordant'
Specify when you use all discordant reads from BAM or CRAM file.
By default, only unmapped reads for mapping to viruses when BAM or CRAM file is specified as an input. If you want to use all discordant reads (e.g. reads without sam flag '2') for mapping to viruses, you can specify '-all_discordant' option. The discordant reads includes read-pairs with distant mapped positions and ones with low MAPQ; therefore, the number of discordant reads are far higher than unmapped reads. This option is only available when '-alignmentin' option was specified. Although this option is less accurate than default setting (= use only unmapped reads), the coverage of reconstructed sequence is longer than default. Please use this at your own risk.

### '-fastqin'
Specify when you use fastq files as inputs. You also need to specify your input file with '-fq1' or '-fq2' options.

### '-fq1 [read_1 fastq file] -fq2 [read_2 fastq file]'
Use when specifing paired fastq files. Available only when '-fastqin' is also specified.

### '-single -fq [fastq file]'
Use when specifying single-end fastq file. Available only when '-fastqin' is also specified.

### '-vref [reference virus genome file]'
Use when specifing reference virus genome file available from NCBI. This option is always required. See 'preparing virus genome reference file' section for detail.

### '-vrefindex [index of reference virus genome file]'
Use when specifing reference virus genome index. This option is always required. See 'preparing virus genome reference file' section for detail.

### '-depth'
When you are using WGS data, you can specify autosome depth with this option. When specified, this pipeline output (average_depth_of_mapped_region / autosome depth). This is particularly important when juding a detected virus sequence is inherited chromosomally integrated sequence or not.

### '-phage'
When specified, this pipeline also reconstructs phage sequences. Generally, there are lots of phage-derived sequences (incld. spike-in PhiX). By default, this pipeline does not reconstruct when a virus genome name contained the word 'phage' or 'Phage'.

### '-picard [picard.jar file]'
Use when specifing the path to 'picard.jar' file. This option is always required. 

### '-p [integer]'
Number of threads. 3 or more is recommended. Default = 2.

### '-bwa'
When you specify this, this tool will use BWA MEM instead of hisat2 when mapping read to virus genomes. You also need to specify bwa index with '-vrefindex' option.

### '-outdir [out_dir_name]'
You can specify arbitrary name for output directory.

### '-overwrite'
In default, this tool does not overwrite directories and files when output directories and files already exist. When you specify this option, this pipeline overwrite existing directories and files.

### '-keep'
You can keep intermediate files when specifying this option.

### '-v', '--version'
Print version.

### '-h', '--help'
Print help message.



# 5. preparing virus genome reference file

###  prepare virus genome reference file
This tool requires a virus reference genome. This file can be downloaded from NCBI.
Here is the instruction to prepare the virus reference genome file.

```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
zcat viral.*.1.genomic.fna.gz > viral_genomic_seq.fa

# Validate the 'viral_genomic_seq.fa' file.
# This file should output two header lines (below). Otherwise, you cannot use this file for this analysis.
# >NC_000898.1 Human herpesvirus 6B, complete genome
# >NC_001664.4 Human betaherpesvirus 6A, variant A DNA, complete virion genome, isolate U1102
cat viral_genomic_seq.fa | grep -e NC_001664.4 -e NC_000898.1
```
Please specify 'viral_genomic_seq.fa' with '-vref' option.

### prepare hisat2 index
When you use hisat2 for mapping (default), you need to make hisat2 index of the virus genome reference file.

```
mkdir hisat2_index
hisat2-build -p [num_thread] viral_genomic_seq.fa ./hisat2_index/viral_genomic_seq
```
Please specify './hisat2_index/viral_genomic_seq' with '-vrefindex' option.

### prepare bwa index
When you use BWA MEM for mapping ('-bwa' option), you need to make bwa index of the virus genome reference file.

```
mkdir bwa_index
bwa index -p ./bwa_index/viral_genomic_seq -a bwtsw viral_genomic_seq.fa
```
Please specify './bwa_index/viral_genomic_seq' with '-vrefindex' option.



# Option. Prepare read depth of autosomes when using a WGS sample.
Before using this script, you need to run 'samtools coverage' to calculate depth of each chromosome. You need to use samtools 1.10 or later to use 'coverage' function. The output from 'samtools coverage' contains mean depth of each chromosome. This program takes the output file from 'samtools coverage' to calculate mean depth of autosomes. You can specify the output file with '-i' option. You can specify names of autosomes by specifying a file containing names of autosomes with -chr option. When you did not specify a file with '-chr' option, this script will use '/path/to/prog/lib/human_autosomes_ucsc_style.txt' by default.
```
samtools coverage --reference ref.fa in.cram > samtools_coverage.txt
python calc_mapping_depth.py -i samtools_coverage.txt -chr /path/to/prog/lib/human_autosomes_ucsc_style.txt -o autosome_depth.txt
```
The output file 'autosome_depth.txt' contains one line with 3 columns. 
- 1st column: average depth of autosomes
- 2nd column: Absolute path of 'samtools_coverage.txt'
- 3rd column: names of the autosomes found in 'samtools_coverage.txt'

You can specify the number in the 1st column with -depth option of the main program.
```
autosome_depth=`cat autosome_depth.txt | cut -f 1`

python main.py \
-depth ${autosome_depth} \
-alignmentin \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
