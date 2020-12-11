#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,pysam
import utils
import log,traceback


def retrieve_unmapped_reads(args, params, filenames):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        # retrieve discordant reads, default
        if args.use_mate_mapped is False and args.all_discordant is False:
            if not args.b is None:
                pysam.view('-@', '%d' % thread_n, '-f', '12', '-F', '3842', '-b', '-o', filenames.discordant_bam, args.b, catch_stdout=False)
            elif not args.c is None:
                pysam.view('-@', '%d' % thread_n, '-f', '12', '-F', '3842', '-b', '-o', filenames.discordant_bam, '--reference', args.fa, args.c, catch_stdout=False)
            pysam.fastq('-@', '%d' % thread_n, '-N', '-0', '/dev/null', '-1', filenames.unmapped_merged_pre1, '-2', filenames.unmapped_merged_pre2, '-s', '/dev/null', filenames.discordant_bam)
            if args.keep is False:
                os.remove(filenames.discordant_bam)
        # retrieve discordant reads, non-default
        else:
            if not args.b is None:
                pysam.view('-@', '%d' % thread_n, '-f', '1', '-F', '3842', '-b', '-o', filenames.discordant_bam, args.b, catch_stdout=False)
            elif not args.c is None:
                pysam.view('-@', '%d' % thread_n, '-f', '1', '-F', '3842', '-b', '-o', filenames.discordant_bam, '--reference', args.fa, args.c, catch_stdout=False)
            pysam.sort('-@', '%d' % thread_n, '-n', '-O', 'BAM', '-o', filenames.discordant_sort_bam, filenames.discordant_bam)
            if args.keep is False:
                os.remove(filenames.discordant_bam)
            if args.all_discordant is True:
                pysam.fastq('-@', '%d' % thread_n, '-N', '-0', '/dev/null', '-1', filenames.unmapped_merged_pre1, '-2', filenames.unmapped_merged_pre2, '-s', '/dev/null', filenames.discordant_sort_bam)
            else:
                pysam.fastq('-@', '%d' % thread_n, '-f', '12', '-F', '3328', '-N', '-0', '/dev/null', '-1', filenames.unmapped_1, '-2', filenames.unmapped_2, '-s', '/dev/null', filenames.discordant_sort_bam)
                if args.use_mate_mapped is True:
                    pysam.view('-@', '%d' % thread_n, '-f', '8', '-F', '3332', '-b', '-o', filenames.unmapped_bam_3, filenames.discordant_sort_bam, catch_stdout=False)
                    pysam.view('-@', '%d' % thread_n, '-f', '4', '-F', '3336', '-b', '-o', filenames.unmapped_bam_4, filenames.discordant_sort_bam, catch_stdout=False)
                    pysam.merge('-@', '%d' % thread_n, '-f', filenames.unmapped_bam_34, filenames.unmapped_bam_3, filenames.unmapped_bam_4)
                    pysam.sort('-@', '%d' % thread_n, '-n', '-O', 'BAM', '-o', filenames.unmapped_sorted_34, filenames.unmapped_bam_34)
                    pysam.fastq('-@', '%d' % thread_n, '-N', '-0', '/dev/null', '-1', filenames.unmapped_3, '-2', filenames.unmapped_4, '-s', '/dev/null', filenames.unmapped_sorted_34)
                # concatenate fastq
                with open(filenames.unmapped_merged_pre1, 'w') as outfile:
                    for f in [filenames.unmapped_1, filenames.unmapped_3]:
                        if os.path.exists(f) is True:
                            with open(f) as infile:
                                for line in infile:
                                    outfile.write(line)
                            utils.gzip_or_del(args, params, f)
                with open(filenames.unmapped_merged_pre2, 'w') as outfile:
                    for f in [filenames.unmapped_2, filenames.unmapped_4]:
                        if os.path.exists(f) is True:
                            with open(f) as infile:
                                for line in infile:
                                    outfile.write(line)
                            utils.gzip_or_del(args, params, f)
        # remove short reads
        infile1=open(filenames.unmapped_merged_pre1)
        infile2=open(filenames.unmapped_merged_pre2)
        outfile1=open(filenames.unmapped_merged_1, 'w')
        outfile2=open(filenames.unmapped_merged_2, 'w')
        min_seq_len=params.min_seq_len
        tmp1,tmp2=[],[]
        for line1,line2 in zip(infile1, infile2):
            tmp1.append(line1)
            tmp2.append(line2)
            if len(tmp1) == 4:
                seqlen1=len(tmp1[1].strip())
                seqlen2=len(tmp2[1].strip())
                if seqlen1 >= min_seq_len and seqlen2 >= min_seq_len:
                    outfile1.write(''.join(tmp1))
                    outfile2.write(''.join(tmp2))
                tmp1,tmp2=[],[]
        infile1.close()
        infile2.close()
        outfile1.close()
        outfile2.close()
        utils.gzip_or_del(args, params, filenames.unmapped_merged_pre1)
        utils.gzip_or_del(args, params, filenames.unmapped_merged_pre2)
        if args.keep is False:
            if os.path.exists(filenames.discordant_sort_bam) is True:
                os.remove(filenames.discordant_sort_bam)
            if args.use_mate_mapped is True:
                os.remove(filenames.unmapped_bam_3)
                os.remove(filenames.unmapped_bam_4)
                os.remove(filenames.unmapped_bam_34)
                os.remove(filenames.unmapped_sorted_34)
    
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

