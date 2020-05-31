#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging


# version
version='2020/05/10'


# HHV-6 refseq IDs
hhv6a_refid='NC_001664.4'
hhv6b_refid='NC_000898.1'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-alignmentin', help='Optional. Specify if you use BAM/CRAM file for input. You also need to specify either -b or -c.', action='store_true')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-all_discordant', help='Optional. Specify if you use all discordant reads from BAM/CRAM file for mapping to viruses. Otherwise, only unmapped reads will be used (default).', action='store_true')
parser.add_argument('-fastqin', help='Optional. Specify if you use unmapped reads for input instead of BAM/CRAM file. You also need to specify -fq1 and -fq2.', action='store_true')
parser.add_argument('-single', help='Optional. Specify if you use single-end unmapped reads for input instead of BAM/CRAM file. Only works when specifing -fastqin option. You also need to specify -fq1.', action='store_true')
parser.add_argument('-fq1', metavar='str', type=str, help='Specify unmapped fastq file, read-1 of read pairs.')
parser.add_argument('-fq2', metavar='str', type=str, help='Specify unmapped fastq file, read-2 of read pairs.')
parser.add_argument('-vref', metavar='str', type=str, help='Required. Specify reference of virus genomes, including HHV-6A and B. Example: viral_genomic_200405.fa')
parser.add_argument('-vrefindex', metavar='str', type=str, help='Required. Specify hisat2 index of virus genomes, including HHV-6A and B. Example: viral_genomic_200405')
parser.add_argument('-depth', metavar='int', type=int, help='Optional. Average depth of input BAM/CRAM file. Only available when using WGS data. If this option is true, will output virus_read_depth/chromosome_depth as well.')
parser.add_argument('-bwa', help='Optional. Specify if you use BWA for mapping instead of hisat2.', action='store_true')
parser.add_argument('-denovo', help='Optional. Specify if you want to perform de-novo assembly.', action='store_true')
parser.add_argument('-picard', metavar='str', type=str, help='Required. Specify full path to picard.jar. Example: /path/to/picard/picard.jar')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 2', default=2)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()


# start
import init
init.init(args, version)


# logging
import log
args.logfilename='for_debug.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
print()
log.logger.info('You are using version "%s"' % version)
log.logger.info('Initial check started.')
initial_check.check(args, sys.argv)


# set up
import setup
setup.setup(args, init.base)
params=setup.params


# output file names
import utils
filenames=utils.empclass()

filenames.discordant_bam      =os.path.join(args.outdir, 'discordant.bam')
filenames.discordant_sort_bam =os.path.join(args.outdir, 'discordant_sorted.bam')
filenames.unmapped_1          =os.path.join(args.outdir, 'unmapped_1.fq')
filenames.unmapped_2          =os.path.join(args.outdir, 'unmapped_2.fq')
filenames.unmapped_3          =os.path.join(args.outdir, 'unmapped_3.fq')
filenames.unmapped_4          =os.path.join(args.outdir, 'unmapped_4.fq')
filenames.unmapped_5          =os.path.join(args.outdir, 'unmapped_5.fq')
filenames.unmapped_6          =os.path.join(args.outdir, 'unmapped_6.fq')
filenames.unmapped_merged_1   =os.path.join(args.outdir, 'unmapped_merged_1.fq')
filenames.unmapped_merged_2   =os.path.join(args.outdir, 'unmapped_merged_2.fq')
filenames.mapped_unsorted_bam =os.path.join(args.outdir, 'mapped_to_virus_orig.bam')
filenames.mapped_sorted       =os.path.join(args.outdir, 'mapped_to_virus_sorted.bam')
filenames.mapped_to_virus_bam =os.path.join(args.outdir, 'mapped_to_virus_dedup.bam')
filenames.mapped_to_virus_bai =os.path.join(args.outdir, 'mapped_to_virus_dedup.bai')
filenames.markdup_metrix      =os.path.join(args.outdir, 'mark_duplicate_metrix.txt')
filenames.bedgraph            =os.path.join(args.outdir, 'mapped_to_virus.bedgraph')
filenames.summary             =os.path.join(args.outdir, 'virus_detection_summary.txt')
filenames.high_cov_pdf        =os.path.join(args.outdir, 'high_coverage_viruses.pdf')

filenames.tmp_bam             =os.path.join(args.outdir, 'tmp.bam')
filenames.tmp_sorted_bam      =os.path.join(args.outdir, 'tmp_sorted.bam')
filenames.tmp_bam_fq1         =os.path.join(args.outdir, 'tmp_bam_1.fq')
filenames.tmp_bam_fq2         =os.path.join(args.outdir, 'tmp_bam_2.fq')
filenames.tmp_rg_bam          =os.path.join(args.outdir, 'tmp_rg.bam')
filenames.tmp_fa              =os.path.join(args.outdir, 'tmp.fa')
filenames.tmp_masked_fa       =os.path.join(args.outdir, 'tmp_masked.fa')
filenames.tmp_fa_dict         =os.path.join(args.outdir, 'tmp.dict')
filenames.hhv6a_vcf_gz        =os.path.join(args.outdir, 'hhv6a.vcf.gz')
filenames.hhv6a_norm_vcf_gz   =os.path.join(args.outdir, 'hhv6a_norm.vcf.gz')
filenames.hhv6a_gatk_naive    =os.path.join(args.outdir, 'hhv6a_reconstructed.fa')
filenames.hhv6b_vcf_gz        =os.path.join(args.outdir, 'hhv6b.vcf.gz')
filenames.hhv6b_norm_vcf_gz   =os.path.join(args.outdir, 'hhv6b_norm.vcf.gz')
filenames.hhv6b_gatk_naive    =os.path.join(args.outdir, 'hhv6b_reconstructed.fa')
filenames.hhv6a_metaspades_o  =os.path.join(args.outdir, 'hhv6a_metaspades_assembly')
filenames.hhv6b_metaspades_o  =os.path.join(args.outdir, 'hhv6b_metaspades_assembly')

# 0. Unmapped read retrieval
if args.alignmentin is True:
    import retrieve_unmapped
    log.logger.info('Unmapped read retrieval started.')
    retrieve_unmapped.retrieve_unmapped_reads(args, params, filenames)
elif args.fastqin is True:
    log.logger.info('Unmapped read retrieval skipped. Read1=%s, read2=%s.' % (args.fq1, args.fq2))
    if args.single is False:
        filenames.unmapped_merged_1=args.fq1
        filenames.unmapped_merged_2=args.fq2
    else:
        filenames.unmapped_merged_1=args.fq1

# 1. mapping
import mapping
log.logger.info('Mapping of unmapped reads started.')
mapping.map_to_viruses(args, params, filenames)
if args.alignmentin is True:
    utils.gzip_or_del(args, params, filenames.unmapped_merged_1)
    utils.gzip_or_del(args, params, filenames.unmapped_merged_2)
log.logger.info('BAM to bedgraph conversion started.')
mapping.bam_to_bedgraph(args, params, filenames)

# 2. identify high coverage viruses
import identify_high_cov
log.logger.info('Identification of high-coverage viruses started.')
identify_high_cov.identify_high_cov_virus_from_bedgraph(args, params, filenames)

# 3. reconstruct HHV-6
import reconstruct_hhv6
if identify_high_cov.hhv6a_highcov is True:
    log.logger.info('HHV-6A sequence reconstruction started.')
    reconstruct_hhv6.reconst_a(args, params, filenames, hhv6a_refid)
if identify_high_cov.hhv6b_highcov is True:
    log.logger.info('HHV-6B sequence reconstruction started.')
    reconstruct_hhv6.reconst_b(args, params, filenames, hhv6b_refid)
if args.keep is False:
    os.remove(filenames.mapped_to_virus_bai)

log.logger.info('All analysis finished!')
