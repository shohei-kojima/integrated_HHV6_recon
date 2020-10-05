#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,subprocess
import log,traceback
import pysam
import utils


def map_to_dr(args, params, filenames, hhv6_refid):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        pysam.view('-bh', '-o', filenames.tmp_bam, filenames.mapped_to_virus_bam, hhv6_refid, catch_stdout=False)
        pysam.sort('-n', filenames.tmp_bam, '-o', filenames.tmp_sorted_bam)
        pysam.fastq('-N', '-0', '/dev/null', '-1', filenames.unmapped_merged_1, '-2', filenames.unmapped_merged_2, '-s', '/dev/null', filenames.tmp_sorted_bam)
        if args.fastqin is True and args.single is True:
            cmd='hisat2 --mp %s -t -x %s -p %d -U %s --no-spliced-alignment | samtools view -Sbh -o %s -' % (params.hisat2_mismatch_penalties, filenames.hhv6_dr_index, thread_n, filenames.unmapped_merged_1, filenames.mapped_unsorted_bam)
        else:
            cmd='hisat2 --mp %s -t -x %s -p %d -1 %s -2 %s --no-spliced-alignment | samtools view -Sbh -o %s -' % (params.hisat2_mismatch_penalties, filenames.hhv6_dr_index, thread_n, filenames.unmapped_merged_1, filenames.unmapped_merged_2, filenames.mapped_unsorted_bam)
        log.logger.debug('mapping command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during mapping.')
            exit(1)
        if not args.keep is True:
            os.remove(filenames.unmapped_merged_1)
            os.remove(filenames.unmapped_merged_2)
        # sort
        pysam.sort('-@', str(thread_n), '-o', filenames.mapped_sorted, filenames.mapped_unsorted_bam)
        if not args.keep is True:
            os.remove(filenames.mapped_unsorted_bam)
        # mark duplicate
        cmd='java -Xms896m -Xmx5376m -jar %s MarkDuplicates CREATE_INDEX=true I=%s O=%s M=%s' % (args.picard, filenames.mapped_sorted, filenames.mapped_to_dr_bam, filenames.markdup_metrix_dr)
        log.logger.debug('picard command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('\n'+ traceback.format_exc())
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        # remove unnecessary files
        os.remove(filenames.tmp_sorted_bam)
        if args.keep is False:
            os.remove(filenames.mapped_sorted)
        # check mapped = 0
        global read_mapped
        read_mapped=True
        with open(filenames.markdup_metrix_dr) as infile:
            for line in infile:
                if 'Unknown Library' in line:
                    ls=line.split()
                    if int(ls[2]) == 0:
                        read_mapped=False
                    break
        
        # convert to bedgraph
        cmd='bamCoverage --outFileFormat bedgraph -p %d --binSize %d -b %s -o %s' % (thread_n, params.bedgraph_bin, filenames.mapped_to_dr_bam, filenames.bedgraph_dr)
        log.logger.debug('bamCoverage command = "'+ cmd +'"')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('\n'+ traceback.format_exc())
            log.logger.error('Error occurred during bamCoverage.')
            exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def output_summary(args, params, filenames):
    log.logger.debug('started.')
    try:
        # load virus names from virus reference seq file
        virus_names={}
        with open(filenames.hhv6_dr_ref) as infile:
            for line in infile:
                if '>' in line:
                    ls=line.strip().split(' ', 1)
                    virus_names[ls[0].replace('>', '')]=ls[1]
        # identify high cov viruses
        prev_id='any'
        high_cov=[]
        for_plot_d={}
        for_plot_cov=[]
        tmp_retain=[]
        with open(filenames.summary_dr, 'w') as outfile:
            with open(filenames.bedgraph_dr) as infile:
                for line in infile:
                    ls=line.split()
                    if ls[0] == prev_id:
                        if int(float(ls[3])) >= 1:
                            cov=int(float(ls[3]))
                            for _ in range(int(ls[1]), int(ls[2])):
                                covs.append(cov)
                    else:
                        if not prev_id == 'any':
                            cov_len=len(covs)
                            genome_covered= cov_len / total_len
                            ave_depth= sum(covs) / total_len
                            if cov_len >= 1:
                                ave_depth_norm= sum(covs) / cov_len
                                if args.depth is not None:
                                    ratio_ave_virus_depth_to_autosome_depth= str(ave_depth_norm / args.depth)
                                else:
                                    ratio_ave_virus_depth_to_autosome_depth='NA'
                                high_cov_judge='False'
                                if genome_covered >= params.genome_cov_thresholds:
                                    if ave_depth_norm >= params.ave_depth_of_mapped_region_threshold:
                                        high_cov.append([prev_id, genome_covered, ave_depth_norm])
                                        for_plot_d[prev_id]=tmp_retain
                                        for_plot_cov.append([ave_depth, prev_id])
                                        high_cov_judge='True'
                            else:
                                ave_depth_norm=0
                                if args.depth is not None:
                                    ratio_ave_virus_depth_to_autosome_depth='0'
                                else:
                                    ratio_ave_virus_depth_to_autosome_depth='NA'
                                high_cov_judge='False'
                            outfile.write('%s_DR\tgenome_length=%d;mapped_length=%d;perc_genome_mapped=%f;average_depth=%f;average_depth_of_mapped_region=%f;ratio_ave_virus_depth_to_autosome_depth=%s\tfasta_header=%s\n' % (prev_id, total_len, cov_len, 100 * genome_covered, ave_depth, ave_depth_norm, ratio_ave_virus_depth_to_autosome_depth, virus_names[prev_id]))
                            tmp_retain=[]
                        total_len=0
                        covs=[]
                        if int(float(ls[3])) >= 1:
                            cov=int(float(ls[3]))
                            for _ in range(int(ls[1]), int(ls[2])):
                                covs.append(cov)
                    total_len += int(ls[2]) - int(ls[1])
                    prev_id=ls[0]
                    tmp_retain.append([ int(i) for i in ls[1:4] ])
            cov_len=len(covs)
            genome_covered= cov_len / total_len
            ave_depth= sum(covs) / total_len
            if cov_len >= 1:
                ave_depth_norm= sum(covs) / cov_len
                if args.depth is not None:
                    ratio_ave_virus_depth_to_autosome_depth= str(ave_depth_norm / args.depth)
                else:
                    ratio_ave_virus_depth_to_autosome_depth='NA'
                high_cov_judge='False'
                if genome_covered >= params.genome_cov_thresholds:
                    if ave_depth_norm >= params.ave_depth_of_mapped_region_threshold:
                        high_cov.append([prev_id, genome_covered, ave_depth_norm])
                        for_plot_d[prev_id]=tmp_retain
                        for_plot_cov.append([ave_depth, prev_id])
                        high_cov_judge='True'
            else:
                ave_depth_norm=0
                if args.depth is not None:
                    ratio_ave_virus_depth_to_autosome_depth='0'
                else:
                    ratio_ave_virus_depth_to_autosome_depth='NA'
                high_cov_judge='False'
            outfile.write('%s_DR\tgenome_length=%d;mapped_length=%d;perc_genome_mapped=%f;average_depth=%f;average_depth_of_mapped_region=%f;ratio_ave_virus_depth_to_autosome_depth=%s\tfasta_header=%s\n' % (prev_id, total_len, cov_len, 100 * genome_covered, ave_depth, ave_depth_norm, ratio_ave_virus_depth_to_autosome_depth, virus_names[prev_id]))
        if len(high_cov) >= 1:
            log.logger.info('high_cov_DR=%s' % ';'.join([ l[0] for l in high_cov ]))
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def mask_low_depth(args, params, filenames, orig_seq_file, refseqid):
    log.logger.debug('started.')
    try:
        depth=[]
        with open(filenames.bedgraph_dr) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] == refseqid:
                    val=int(ls[3])
                    for _ in range(int(ls[1]), int(ls[2])):
                        depth.append(val)
        orig_seq=[]
        with open(orig_seq_file) as infile:
            for line in infile:
                if '>' in line:
                    header=line.strip()
                    header += ' masked\n'
                else:
                    orig_seq.append(line.strip())
        orig_seq=''.join(orig_seq)
        if not len(depth) == len(orig_seq):
            log.logger.error('Error occurred during making low depth seq. Length of refseq and depth info is not equal.')
            exit(1)
        masked_seq=[]
        for seq,val in zip(orig_seq, depth):
            if val >= params.reconst_minimum_depth:
                masked_seq.append(seq)
            else:
                masked_seq.append('N')
        masked_seq.append('\n')
        with open(filenames.tmp_masked_fa, 'w') as outfile:
            outfile.write(header)
            outfile.write(''.join(masked_seq))
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)



def reconst_a(args, params, filenames, refseqid):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        if args.alignmentin is True:
            sample_name=os.path.basename(args.b) if not args.b is None else os.path.basename(args.c)
        else:
            sample_name=os.path.basename(args.fq1)
        pysam.view(filenames.mapped_to_dr_bam, '-h', '-o', filenames.tmp_bam, refseqid, catch_stdout=False)
        _,seq=utils.retrieve_only_one_virus_fasta(filenames.hhv6_dr_ref, refseqid)
        with open(filenames.tmp_fa, 'w') as outfile:
            outfile.write('>%s DR %s\n%s\n' % (refseqid, sample_name, seq))
        pysam.faidx(filenames.tmp_fa)
        
        # mask low depth regions
        mask_low_depth(args, params, filenames, filenames.tmp_fa, refseqid)
        
        if os.path.exists(filenames.tmp_fa_dict) is True:
            os.remove(filenames.tmp_fa_dict)
        cmd='java -jar %s CreateSequenceDictionary R=%s O=%s' % (args.picard, filenames.tmp_fa, filenames.tmp_fa_dict)
        log.logger.debug('picard command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        cmd='java -jar %s AddOrReplaceReadGroups I=%s O=%s RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20' % (args.picard, filenames.tmp_bam, filenames.tmp_rg_bam)
        log.logger.debug('picard command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
#        pysam.index('-@', str(thread_n), filenames.tmp_rg_bam)
        pysam.index(filenames.tmp_rg_bam)
        cmd='gatk --java-options "-Xmx4g" HaplotypeCaller -R %s -I %s -O %s' % (filenames.tmp_fa, filenames.tmp_rg_bam, filenames.hhv6a_dr_vcf_gz)
        log.logger.debug('gatk command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        cmd='bcftools norm -c x -f %s %s -Oz -o %s' % (filenames.tmp_masked_fa, filenames.hhv6a_dr_vcf_gz, filenames.hhv6a_dr_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools index %s' % filenames.hhv6a_dr_norm_vcf_gz
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools consensus -f %s -o %s %s' % (filenames.tmp_masked_fa, filenames.hhv6a_dr_gatk_naive, filenames.hhv6a_dr_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        # remove unnecessary files
        os.remove(filenames.tmp_rg_bam)
        os.remove(filenames.tmp_rg_bam +'.bai')
        os.remove(filenames.tmp_fa)
        os.remove(filenames.tmp_fa +'.fai')
        os.remove(filenames.tmp_masked_fa)
        os.remove(filenames.tmp_masked_fa +'.fai')
        os.remove(filenames.tmp_fa_dict)
        if args.keep is False:
            os.remove(filenames.hhv6a_dr_vcf_gz +'.tbi')
            os.remove(filenames.hhv6a_dr_norm_vcf_gz)
            os.remove(filenames.hhv6a_dr_norm_vcf_gz +'.csi')
        # remove unnecessary files
        os.remove(filenames.tmp_bam)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

def reconst_b(args, params, filenames, refseqid):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        if args.alignmentin is True:
            sample_name=os.path.basename(args.b) if not args.b is None else os.path.basename(args.c)
        else:
            sample_name=os.path.basename(args.fq1)
        pysam.view(filenames.mapped_to_dr_bam, '-h', '-o', filenames.tmp_bam, refseqid, catch_stdout=False)
        _,seq=utils.retrieve_only_one_virus_fasta(filenames.hhv6_dr_ref, refseqid)
        with open(filenames.tmp_fa, 'w') as outfile:
            outfile.write('>%s DR %s\n%s\n' % (refseqid, sample_name, seq))
        pysam.faidx(filenames.tmp_fa)
        
        # mask low depth regions
        mask_low_depth(args, params, filenames, filenames.tmp_fa, refseqid)
        
        if os.path.exists(filenames.tmp_fa_dict) is True:
            os.remove(filenames.tmp_fa_dict)
        cmd='java -jar %s CreateSequenceDictionary R=%s O=%s' % (args.picard, filenames.tmp_fa, filenames.tmp_fa_dict)
        log.logger.debug('picard command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        cmd='java -jar %s AddOrReplaceReadGroups I=%s O=%s RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20' % (args.picard, filenames.tmp_bam, filenames.tmp_rg_bam)
        log.logger.debug('picard command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
#        pysam.index('-@', str(thread_n), filenames.tmp_rg_bam)
        pysam.index(filenames.tmp_rg_bam)
        cmd='gatk --java-options "-Xmx4g" HaplotypeCaller -R %s -I %s -O %s' % (filenames.tmp_fa, filenames.tmp_rg_bam, filenames.hhv6b_dr_vcf_gz)
        log.logger.debug('gatk command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        cmd='bcftools norm -c x -f %s %s -Oz -o %s' % (filenames.tmp_masked_fa, filenames.hhv6b_dr_vcf_gz, filenames.hhv6b_dr_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools index %s' % filenames.hhv6b_dr_norm_vcf_gz
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools consensus -f %s -o %s %s' % (filenames.tmp_masked_fa, filenames.hhv6b_dr_gatk_naive, filenames.hhv6b_dr_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        # remove tmp files
        os.remove(filenames.tmp_rg_bam)
        os.remove(filenames.tmp_rg_bam +'.bai')
        os.remove(filenames.tmp_fa)
        os.remove(filenames.tmp_fa +'.fai')
        os.remove(filenames.tmp_masked_fa)
        os.remove(filenames.tmp_masked_fa +'.fai')
        os.remove(filenames.tmp_fa_dict)
        # remove unnecessary files
        if args.keep is False:
            os.remove(filenames.hhv6b_dr_vcf_gz +'.tbi')
            os.remove(filenames.hhv6b_dr_norm_vcf_gz)
            os.remove(filenames.hhv6b_dr_norm_vcf_gz +'.csi')
        # remove unnecessary files
        os.remove(filenames.tmp_bam)
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

