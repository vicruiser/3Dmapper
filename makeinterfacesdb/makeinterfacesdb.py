# coding: utf-8

# Import necesary modules
import sys
import os
import gzip
import re
import emoji
import glob
import subprocess
import vcfpy
import time
import os.path

import pandas as pd
import numpy as np

from halo import Halo
from timeit import default_timer as timer
from subprocess import call

# import functions from scripts
from .parse_argv import parse_commandline
from .split import split
from .decorator import tags
from .logger import get_logger


def main():
    # aesthetics
    description = '''

    ----------------------------------------- Welcome to ----------------------------------------------

    $$$$$$$\  $$$$$$$\  $$$$$$$\   
    $$  __$$\ $$  __$$\ $$  __$$\     
    $$ |  $$ |$$ |  $$ |$$ |  $$ |$$$$$$\$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\ 
    $$$$$$$  |$$ |  $$ |$$$$$$$\ |$$  _$$  _$$\  \____$$\ $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$\ 
    $$  ____/ $$ |  $$ |$$  __$$\ $$ / $$ / $$ | $$$$$$$ |$$ /  $$ |$$ /  $$ |$$$$$$$$ |$$ |  \__|
    $$ |      $$ |  $$ |$$ |  $$ |$$ | $$ | $$ |$$  __$$ |$$ |  $$ |$$ |  $$ |$$   ____|$$ |
    $$ |      $$$$$$$  |$$$$$$$  |$$ | $$ | $$ |\$$$$$$$ |$$$$$$$  |$$$$$$$  |\$$$$$$$\ $$ |
    \__|      \_______/ \_______/ \__| \__| \__| \_______|$$  ____/ $$  ____/  \_______|\__|
                                                            $$ |      $$ |
                                                            $$ |      $$ |
                                                            \__|      \__|

    ---------------  Map annotated genomic variants to protein interfaces data in 3D. -----------------

    \n'''

    epilog = \
        '''
          -------------------------------------------------------------------------        
         |  Copyright (c) 2019 Victoria Ruiz --                                    |  
         |  vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          -------------------------------------------------------------------------

        '''
    # print ascii art
    print(description)
    print(epilog)
    # initialize spinner decorator
    spinner = Halo(text='Loading', spinner='dots12', color="red")
    # parse command line options
    args = parse_commandline()

    # set out dir and out file names
    # created by default
    out_dir = os.path.join(args.out, 'DBs')
    # out_file = os.path.join(
    #    out_dir, 'variants.vep')  # created by default
    # set output dir to split vep
    intdb_outdir = os.path.join(out_dir, 'intDB')  # created by default
    # create output dir if it doesn't exist
    os.makedirs(intdb_outdir, exist_ok=True)

    # set up the report
    report = open(os.path.join(out_dir, 'makeinterfacesdb.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('''
    Command line input: 
    -------------------
    \n''')
    report.write((" ".join(sys.argv)) + '\n' + '\n' + '\n')
    time_format = '[' + time.ctime(time.time()) + '] '
    start = time.time()

    if args.intdb is not None:
        # set up a log file
        logger = get_logger('main', out_dir)
        log_dir = out_dir
        logger.info('Reading and splitting input file.')

        # report info
        report.write(time_format + 'Reading and splitting input file. \n')
        # for loop in case we have multiple inputs to read from a list of files
        for f in args.intdb:
            # check if input is a file
            if isfile(f) == 'list_files':
                with open(f) as list_int_files:
                    int_f = list_int_files.read().splitlines()
                    logger.info(
                        'Input interfaces file contains a list of variants files to process.')
                    # for every prot id
                    for int_infile in int_f:
                        # split interface db
                        split('ENSP', int_infile, intdb_outdir,
                              'txt', args.force, log_dir)
                        # log info
                        logger.info(
                            int_infile + ' has been splitted successfully.')

            elif isfile(f) == 'is_file':
                # split interface db
                split('ENSP', f, intdb_outdir,
                      'txt', args.force, log_dir)
                # log info
                logger.info(
                    f + ' has been splitted successfully.')
            else:
                logger.error(
                    'Error: Interfaces file input could not be found.')
                raise IOError

            # finish report
            end = time.time()
            report.write(
                time_format + 'Reading and splitting input file...done. \n')
            report.write(
                time_format + 'Generation of interfaces DB in ' + intdb_outdir + ' took ' + str(round(end-start, 2)) + 's\n')
            report.close()
