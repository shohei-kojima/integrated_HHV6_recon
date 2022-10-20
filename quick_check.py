#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging


# version
version='2022/10/20'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')
parser.add_argument('-bl', metavar='str', type=str, help='Either -b or -c is Required. Specify file that contains paths to input paired-end BAM files.')
parser.add_argument('-cl', metavar='str', type=str, help='Either -b or -c is Required. Specify file that contains paths to input paired-end CRAM files.')
parser.add_argument('-s', metavar='str', type=str, help='Specify file that contains unmapped SAM files.')
parser.add_argument('-sl', metavar='str', type=str, help='Specify file that contains file paths of unmapped SAM files.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_quick_check', default='./result_quick_check')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 3 or more is recommended. Default: 1', default=1)
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
parser.add_argument('-singularity', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()


# start
import init
init.init(args, version)


# logging
import log
args.logfilename='for_debug_quick_check.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
print()
log.logger.info('You are using version "%s"' % version)
log.logger.info('Initial check started.')
initial_check.check_quick_check(args, sys.argv, init.base)


# set up
import setup
setup.setup(args, init.base)
params=setup.params


# output file names
import utils
filenames=utils.empclass()

filenames.unmapped       =os.path.join(args.outdir, 'unmapped.fq')
filenames.mapped_sam     =os.path.join(args.outdir, 'mapped_to_hhv6.sam')
filenames.final_result   =os.path.join(args.outdir, 'HHV6_check_results.txt')


# quick check
import quick_checking
log.logger.info('Quick checking started.')
quick_checking.load_files(args, params, filenames)
quick_checking.checking(args, params, filenames)

log.logger.info('Quick checking finished!')
