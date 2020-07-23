#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,multiprocessing
from os.path import abspath,dirname,realpath,join
import log,traceback


def setup(args, base):
    log.logger.debug('started')
    try:
        sys.setrecursionlimit(10000)
        
        # load parameter settings
        import load_parameters
        global params
        params=load_parameters.load(args)

    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
