# This is a README file for iciHHV6_reconstruction_pipeline

# Flowchart
<img src='https://github.com/shohei-kojima/iciHHV6_reconstruction/blob/master/lib/image_for_readme.png' width='320px'>


# 1. Getting started
## 1.1 Use Singularity container (recommended)
This tool is available as a Singularity container.
Currently, the Singularity container is provided directly from Genome Immunobiology Lab, RIKEN.
Please contact us when you would like to use this.

It requires Singularity version 3 (does not work on version 2) to run the container.
Please see [here](https://sylabs.io/guides/3.0/user-guide/installation.html) for installation of Singularity.

You do NOT need to set up environment described in the next chapeter `1.2 Use scripts available from GitHub`, when you use the Singularity container.

## 1.2 Use scripts available from GitHub
You can use the original scripts available from GitHub instead of the Singularity container.
In this case, you need to set up environment in your Linux for this tool.
Envriomnet set-up described below is not needed if you use the Singularity container.
We tested this tool in Ubuntu 18.04, CentOS 7.7.1908, Red Hat 7.3.

### Required software for running with default settings
All tools listed here are needed to be installed and added to $PATH.
We recommend to use Anaconda3.

- Python 3.7 or later  **Please use Python 3.7 or later. Python 3.6 does NOT work!**
- hisat2 2.1.0 or later 
- samtools 1.10 or later
- bcftools 1.9 or later
- bamCoverage 3.3.1 or later (deeptools 3.3.1 or later)
- gatk 4.1.7.0 or later
- picard 2.21.9 or later
- java (version matches to your picard)

### Optional prerequisites
- bwa 0.7.17 or later (if you specify '-bwa' option)
- SPAdes genome assembler v3.13.1 or later (if you specify '-denovo' option)

### Required Python modules
- pysam 0.15.2 or later
- matplotlib 3.1.1 or later
- built-in modules: os,sys,datetime,multiprocessing,logging,traceback,argparse,glob,pylab,subprocess,gzip

### Download this tool from GitHub (currently not public)
git clone https://github.com/shohei-kojima/iciHHV6_reconstruction

### Install python modules
```
conda install -c bioconda pysam 
conda install -c bioconda matplotlib
```

### Install hisat2, bcftools, deeptools
```
conda install -c bioconda hisat2 
conda install -c bioconda bcftools 
conda install -c bioconda deeptools
```

### Install samtools
Because this pipeline requires samtools1.10, please do NOT use `conda` for installation of samtools. Conda installs samtools 1.9. Please download a latest version from htslib homepage (http://www.htslib.org/download/). 
```
# install
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xf samtools-1.10.tar.bz2
cd samtools-1.10
./configure --prefix=/where/to/install
make
make install

# add below to '~/.bashrc'
export PATH=/where/to/install/bin:$PATH
```

### Install gatk
```
# install
wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.0/gatk-4.1.8.0.zip
unzip gatk-4.1.8.0.zip

# add below to '~/.bashrc'
export PATH=/where/to/install/gatk-4.1.8.0:$PATH
```

### Install picard
```
# install
wget https://github.com/broadinstitute/picard/releases/download/2.23.1/picard.jar

# test run; command below should print help message if picard works without problem
java -jar /where/to/install/picard.jar -h
```


# 2. quick usage guide for the impatient
## 2.1 Run Singularity container (recommended)
We highly recommend to use Singularity container. This container have exactly the same functionality to the original scripts in GitHub.
This tool takes a BAM or CRAM file containing PAIRED-end reads. BAM/CRAM file with single-end reads is NOT supported.

### when you use your BAM file as an input
```
singularity exec iciHHV6_reconstruction_[version].sif reconst \
-alignmentin \
-b your_file.bam \
-p 4
```

### when you use your CRAM file as an input
```
singularity exec iciHHV6_reconstruction_[version].sif reconst \
-alignmentin \
-c your_file.cram \
-fa /path/to/reference_genome_of_cram.fa \
-p 4
```

### when you use your fastq files as inputs
```
singularity exec iciHHV6_reconstruction_[version].sif reconst \
-fastqin \
-fq1 your_file_1.fastq \
-fq2 your_file_2.fastq \
-p 4
```


## 2.2 Run Script downloaded from GitHub
This tool is also available from GitHub. When Singularity is not available, you can use the scripts from GitHub.
In this case, you need to set up environment required for this tool. For details, please see `1.2 Set up required environment in your Linux`. 

### when you use your BAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-b your_file.bam \
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
-c your_file.cram \
-fa /path/to/reference_genome_of_cram.fa \
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
-b your_file.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```

This tool takes a BAM or CRAM file containing PAIRED-end reads. BAM/CRAM file with single-end reads is NOT supported.

By default, only unmapped reads are used for mapping to viruses when BAM or CRAM file is specified as an input. If you want to use all discordant reads (e.g. reads without sam flag '2') for mapping to viruses, you can specify '-all_discordant' option. Discordant reads include read-pairs with distantly mapped positions and ones with low MAPQ; therefore, the number of discordant reads are far higher than unmapped reads. This option is only available when '-alignmentin' option is specified. This option is suspected to be less accurate than using default settingss (= using only unmapped reads), however the coverage of reconstructed sequence may be higher than default. Please be aware of this caveat when using this option.

### when you use your fastq files as inputs (fastqin option)
```
python main.py \
-fastqin \
-fq1 your_file_1.fastq \
-fq2 your_file_2.fastq \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_hisat2_index \
-picard /path/to/picard.jar \
-p 4
```
You can specify your fastq files as inputs by '-fastqin' option. This tool does NOT support single-end reads. Please specify paired-end reads as separate files ('-fq1' and '-fq2' options).

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
'-denovo' option performs de novo assembly of the HHV-6 sequence, in addition to variant calling and consensus sequence generation.



# 3. output files
### Here are the typical output files
```
# output files
./result_out/
├── for_debug.log
├── hhv6a_reconstructed.fa
├── hhv6a.vcf.gz
├── high_coverage_viruses.pdf
├── mapped_to_virus.bedgraph
├── mapped_to_virus_dedup.bam
├── mark_duplicate_metrix.txt
└── virus_detection_summary.txt

# file size (30x BAM was used as input)
24K     ./result_out/for_debug.log
160K    ./result_out/hhv6a_reconstructed.fa
172K    ./result_out/hhv6a.vcf.gz
400K    ./result_out/high_coverage_viruses.pdf
1.2M    ./result_out/mapped_to_virus.bedgraph
8.2M    ./result_out/mapped_to_virus_dedup.bam
4.0K    ./result_out/mark_duplicate_metrix.txt
3.3M    ./result_out/virus_detection_summary.txt

# brief explanation of output files
./result_out/
├── for_debug.log               # log file
├── hhv6a_reconstructed.fa      # reconstructed HHV-6
├── hhv6a.vcf.gz                # VCF of HHV-6
├── high_coverage_viruses.pdf   # plot of HHV-6 read depth
├── mapped_to_virus.bedgraph    # mapping coverage
├── mapped_to_virus_dedup.bam   # bam with reads used for analysis
├── mark_duplicate_metrix.txt   # gatk markduplicate output
└── virus_detection_summary.txt # main summary file
```


### 'virus_detection_summary.txt'
This is one of the main result files for most users. This contains read coverage information of each virus genome used in the viral genomes input file. Currently, this pipeline is supported for detecting and reconstructing endogenous or chromosomally-integrated HHV-6, in which case exogenous HHV-6A and HHV-6B reference sequences are used. However, the computational cost of adding additional viral genomes for detection is low, and we have been able to detect (but not reconstruct) other integrated viruses using this pipeline.

- 1st column: RefSeq ID (fasta header without attribution)
- 2nd column: Whether virus-derived reads are abundantly detected (consistent with germline or high-fraction mosaic integration) or not
- 3rd column: Coverage information
    - genome_length: Length of the virus genome
    - mapped_length: Length of the virus genome covered by one or more mapped reads
    - perc_genome_mapped: Percent of virus genome with one or more reads
    - average_depth: Average of mapped read depth (average of whole viral genome)
    - average_depth_of_mapped_region: Average of mapped read depth of mapped regions (average of only mapped regions)
    - ratio_ave_virus_depth_to_autosome_depth: Average of mapped read depth of mapped regions divided by autosome depth provided with '-depth' option. Available only when '-depth' option was specified. Otherwise, 'NA'.
- 4th column: attribution of fasta header

### 'hhv6a_reconstructed.fa', 'hhv6b_reconstructed.fa'
This is one of the main result files for most users. This file contains HHV-6 sequence reconstructed with called variants. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B. Genomic regions where do not have any reads (= 0 reads mapped) are masked by a character 'N.'

### 'high_coverage_viruses.pdf'
This is one of the main result files for most users who are using this tool to screen for additional potential integrated viruses. If there are one or more viruses in your sample, this tool outputs read coverage of those viruses for visual inspection.

### 'hhv6a.vcf.gz', 'hhv6b.vcf.gz'
This is one of the main result files for most users. This file contains variant calls. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B.

### 'mapped_to_virus_dedup.bam'
BAM file containing alignment with virus genomes.

### 'mapped_to_virus.bedgraph'
Read depth of virus genomes.

### 'mark_duplicate_metrix.txt'
Summary of picard MarkDuplicates running with 'mapped_to_virus_dedup.bam.'

### 'hhv6a_metaspades_assembly', 'hhv6b_metaspades_assembly'
This file contains the integrated HHV-6 sequence reconstructed by metaspades. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B. You need to specify '-denovo' option to obtain this result.

### 'for_debug.log'
Log file. Stores all information, including input arguments and parameter setting.



# 4. options
### '-alignmentin'
Specify when you use a BAM or CRAM file as an input. You need to specify which input file type you are using with '-b' or '-c' option. When you use CRAM file as an input, you also need to specify reference genome with '-fa' option.

### '-b [BAM file]'
Use when specifing BAM file. Available only when '-alignmentin' is also specified.

### '-c [CRAM file] -fa [reference fasta file]'
Use when specifing CRAM file. Available only when '-alignmentin' is also specified.

### '-all_discordant'
Specify when you want to use all discordant reads from a BAM or CRAM file.
By default, only unmapped reads for are used for mapping to viruses when a BAM or CRAM file is specified as the input. If you want to use all discordant reads (e.g. reads without sam flag '2') for mapping to viruses, you can specify '-all_discordant' option. The discordant reads includes read-pairs with distantly mapped positions and ones with low MAPQ; therefore, the number of discordant reads are far higher than unmapped reads. This option is only available when '-alignmentin' option is specified. Please be aware of potential artifacts due to use of this option, e.g. mapping to repeatitive sequences present in both the human reference and the viral genome.

### '-fastqin'
Specify when you use fastq files as input. You also need to specify your input files with '-fq1' or '-fq2' options.

### '-fq1 [read_1 fastq file] -fq2 [read_2 fastq file]'
Use when specifing paired fastq files. Available only when '-fastqin' is also specified.

### '-single -fq [fastq file]'
Use when specifying single-end fastq file. Available only when '-fastqin' is also specified.

### '-vref [reference virus genome file]'
You do not need to specify this when running the Singularity container. Singularity container already contains iciHHV-6A and iciHHV-6B reference sequences (NC_001664.4 and NC_000898.1, respectively).
Use when specifing reference virus genome file available from NCBI. This option is always required when running scripts downloaded from GitHub. See 'preparing virus genome reference file' section for details.

### '-vrefindex [index of reference virus genome file]'
You do not need to specify this when running the Singularity container. Singularity container already contains hisat2 index of iciHHV-6A and iciHHV-6B reference sequences (NC_001664.4 and NC_000898.1, respectively).
Use when specifing reference virus genome index. This option is always required when running scripts downloaded from GitHub. See 'preparing virus genome reference file' section for detail.

### '-depth'
This is optional, but we highly recommend to use this. When you are using WGS data, you can specify autosome depth with this option. When specified, this pipeline outputs a value (average_depth_of_mapped_region / autosome depth). This is particularly important when juding whether a detected virus sequence is present with most cell's DNA or not. Please see a section [6. (Optional) Prepare read depth of autosomes when using a WGS sample.] for more details.

### '-phage'
When specified, this pipeline also reconstructs phage sequences. Generally, there are lots of phage-derived sequences (including spike-in PhiX) in  WGS data. By default, this pipeline does not reconstruct virus genomes whose name contains the word 'phage' or 'Phage'.

### '-picard [picard.jar file]'
Use when specifing the path to 'picard.jar' file. This option is always required. 

### '-p [integer]'
Number of threads. 3 or more is recommended. Default = 2.

### '-bwa'
When you specify this, this tool will use BWA MEM instead of hisat2 when mapping read to virus genomes. You also need to specify bwa index with '-vrefindex' option.

### '-outdir [out_dir_name]'
You can specify your chosen name for the output directory.

### '-overwrite'
In default, this tool does not overwrite directories and files when output directories and files already exist. When you specify this option, this pipeline  will overwrite existing directories and files.

### '-keep'
You can keep intermediate files when specifying this option.

### '-v', '--version'
Print version.

### '-h', '--help'
Print help message.



# 5. preparing virus genome reference file (required only when running GitHub scripts, not required to run the Singularity container)

###  prepare virus genome reference file
This tool requires a virus reference genome. This file can be downloaded from NCBI.
Here are instructions for how to prepare the virus reference genome file.

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



# 6. (Optional) Prepare read depth of autosomes when using a WGS sample.
We highly recommend to input autosome depth of an input WGS to calculate a value: average_depth_of_mapped_region / autosome depth.
You first need to run `samtools coverage` to calculate depth of each chromosome. The output from `samtools coverage` contains mean depth of each chromosome. This program takes the output file from `samtools coverage` to calculate mean depth of autosomes. Please specify the output file from `samtools coverage` with the '-i' option. Please specify names of autosomes by specifying a file containing names of autosomes with '-chr' option. When you did not specify a file with '-chr' option, this script will use '/path/to/prog/lib/human_autosomes_ucsc_style.txt' by default.
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

# LICENSE
Copyright (c) 2020 RIKEN
All Rights Reserved.
See file LICENSE for details.

This tool was developed by Shohei Kojima and Nicholas F. Parrish.
