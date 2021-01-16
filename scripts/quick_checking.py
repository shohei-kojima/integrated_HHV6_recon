#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''

import os,sys,subprocess,pysam
import utils
import log,traceback


def load_files(args, params, filenames):
    log.logger.debug('started.')
    try:
        def bam_check(path):
            if os.path.exists(path) is False:
                log.logger.error('Input file (%s) was not found.' % path)
                exit(1)
            if os.path.exists(path +'.bai') is False:
                if os.path.exists(path[:-4] +'bai') is False:
                    log.logger.error('Index file for input file (%s) was not found.' % path)
                    exit(1)
                else:
                    index_time=os.stat(path[:-4] +'bai').st_mtime
            else:
                index_time=os.stat(path +'.bai').st_mtime
            diff= os.stat(path).st_mtime - index_time
            if diff > 0:
                log.logger.error('BAM index is older than BAM. Please generate the index again.')
                exit(1)

        def cram_check(path):
            if os.path.exists(path) is False:
                log.logger.error('Input file (%s) was not found.' % path)
                exit(1)
            if os.path.exists(path +'.crai') is False:
                if os.path.exists(path[:-4] +'crai') is False:
                    log.logger.error('Index file for input file (%s) was not found.' % path)
                    exit(1)
                else:
                    index_time=os.stat(path[:-4] +'crai').st_mtime
            else:
                index_time=os.stat(path +'.crai').st_mtime
            diff= os.stat(path).st_mtime - index_time
            if diff > 0:
                log.logger.error('CRAM index is older than CRAM. Please generate the index again.')
                exit(1)
            
        filenames.fpaths=[]
        if args.b is not None:
            bam_check(args.b)
            filenames.fpaths.append(args.b)
            args.file_type='rb'
        elif args.c is not None:
            cram_check(args.c)
            filenames.fpaths.append(args.c)
            args.file_type='rc'
        elif args.bl is not None:
            with open(args.bl) as infile:
                for line in infile:
                    path=line.strip()
                    bam_check(path)
                    filenames.fpaths.append(path)
            args.file_type='rb'
        elif args.cl is not None:
            with open(args.cl) as infile:
                for line in infile:
                    path=line.strip()
                    cram_check(path)
                    filenames.fpaths.append(path)
            args.file_type='rc'
        log.logger.info('%d file(s) will be analyzed.' % len(filenames.fpaths))
        
    except SystemExit:
        log.logger.debug('\n'+ traceback.format_exc())
        exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
    

def checking(args, params, filenames):
    log.logger.debug('started.')
    try:
        read_num_limit=params.quick_check_read_num
        no_hhv_threshold= 3 / 1000000
        need_check_threshold= 20 / 1000000
        dr_threshold= 150 / 1000000
        n_false=0
        n_need_check=0
        n_dr=0
        n_full=0
        
        finalfile=open(filenames.final_result, 'w')
        finalfile.write('#file\tnum_unmapped_read_analyzed\tnum_read_mapped_to_HHV6\tHHV6_exists?\n')
        for f in filenames.fpaths:
            if args.file_type == 'rb':
                infile=pysam.AlignmentFile(f, 'rb', check_sq=False)
            elif args.file_type == 'rc':
                infile=pysam.AlignmentFile(f, 'rc', reference_filename=args.fa)
            n=0
            with open(filenames.unmapped, 'w') as outfile:
                tmp=[]
                for read in infile.fetch('*', until_eof=True):
                    if read.is_unmapped:
                        if not 'TAACCC' in read.query_sequence and not 'GGGTTA' in read.query_sequence:
                            if read.is_read1 is True:
                                header='@%s/1' % read.query_name
                            else:
                                header='@%s/2' % read.query_name
                            tmp.append('%s\n%s\n+\n%s\n' % (header, read.query_sequence, read.qual))
                            n += 1
                    if len(tmp) == 100_000:
                        outfile.write(''.join(tmp))
                        tmp=[]
                    if n == read_num_limit:
                        break
                if len(tmp) >= 1:
                    outfile.write(''.join(tmp))
                outfile.flush()
                os.fdatasync(outfile.fileno())
            infile.close()
            if n == 0:
                log.logger.info('No unmapped reads found in %s. Will continue anyway.' % (n, f))
                finalfile.write('%s\t%d\tNA\tNA\n' % (f, n))
                utils.gzip_or_del(args, params, filenames.unmapped)
                continue
            elif n < read_num_limit:
                log.logger.warning('Only %d unmapped reads were found in %s. Will continue anyway.' % (n, f))
            # mapping
            cmd='hisat2 --mp %s -t -x %s -p %d -U %s --no-spliced-alignment > %s' % (params.hisat2_mismatch_penalties, args.vrefindex, args.p, filenames.unmapped, filenames.mapped_sam)
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during mapping.')
                exit(1)
            utils.gzip_or_del(args, params, filenames.unmapped)
            # count mapped
            mapped_n=0
            with open(filenames.mapped_sam) as infile:
                for line in infile:
                    if not line[0] == '@':
                        ls=line.split()
                        if not ls[5] == '*':
                            readlen=len(ls[9])
                            if ls[5] == '%dM' % readlen:
                                mapped_n += 1
            mapped_ratio= mapped_n / n
            if mapped_ratio < no_hhv_threshold:
                judge='False'
                n_false += 1
            elif mapped_ratio < need_check_threshold:
                judge='Need_further_check'
                n_need_check += 1
            elif mapped_ratio < dr_threshold:
                judge='likely_solo-DR'
                n_dr += 1
            else:
                judge='likely_Full-length'
                n_full += 1
            finalfile.write('%s\t%d\t%d\t%s\n' % (f, n, mapped_n, judge))
        utils.gzip_or_del(args, params, filenames.mapped_sam)
        finalfile.flush()
        os.fdatasync(finalfile.fileno())
        log.logger.info('\n\n\033[34mQuick check result:\n\n  No HHV-6 = %d\n  Need check = %d\n  Likely solo-DR = %d\n  Likely Full-length = %d\033[0m\n\n  \033[31mCaveats: This result is estimation and only for a screening purpose. This is not a conclusive result.\033[0m\n' % (n_false, n_need_check, n_dr, n_full))
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
