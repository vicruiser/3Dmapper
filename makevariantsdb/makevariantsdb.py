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
import logging

import os.path as path
import pandas as pd
import numpy as np

from halo import Halo
from timeit import default_timer as timer
from subprocess import call
from multiprocessing import Pool

# import functions from scripts
from .parse_argv import parse_commandline
from .run_vep import run_vep
from .split import split
from .detect_vcf_format import detect_format
from .vcf2vep import vcf2vep
from .maf2vep import maf2vep
from .add_header import add_header
from .decorator import tags
from .logger import get_logger
from .input_isfile import isfile

# num_cpus = psutil.cpu_count(logical=False)


class generateVarDB:

    def vcf(self, var_infile, out_dir, out_file, vardb_outdir, overwrite, log_dir, parallel=False):

        # from vcf to vep
        vcf2vep(var_infile, out_dir,
                out_file, overwrite, log_dir, parallel)
        # add header to resulting vep file
        add_header(out_file)

        # split vep file by protein id to speed up the
        # mapping process
        split('ENSG', out_file, vardb_outdir,
              'vep', overwrite, log_dir, parallel)

    def vep(self, var_infile, vardb_outdir, overwrite, log_dir, parallel=False):

        # split vep file by protein id to speed up the
        # mapping process
        split('ENSG', var_infile, vardb_outdir,
              'vep', overwrite, log_dir, parallel)

    def maf(self, var_infile, out_dir, out_file, vardb_outdir, overwrite, log_dir, parallel=False):
        # from vcf to vep
        maf2vep(var_infile, out_dir,
                out_file, overwrite, log_dir)
        # split vep file by protein id to speed up the
        # mapping process
        split('ENSG', out_file, vardb_outdir,
              'vep', overwrite, log_dir, parallel)


def main():
    # parse command line options
    args = parse_commandline()
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
    # set out dir and out file names
    # created by default
    out_dir = os.path.join(args.out, 'DBs')
    out_file = os.path.join(
        out_dir, 'variants.vep')  # created by default
    # set output dir to split vep
    vardb_outdir = os.path.join(out_dir, 'varDB')  # created by default
    # create output dir if it doesn't exist
    os.makedirs(vardb_outdir, exist_ok=True)
    # initialize class
    varfile = generateVarDB()

    # set up the results report
    report = open(os.path.join(out_dir, 'makevariantsdb.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    report.write((" ".join(sys.argv)) + '\n' + '\n' + '\n')
    time_format = '[' + time.ctime(time.time()) + '] '

    # set up a log file
    logger = get_logger('main', out_dir)
    log_dir = out_dir
    # spinner.start(text=' Running makevariantsdb...')
    start = time.time()

    # change input format if file doesn't exists or overwrite is True
    if not os.listdir(vardb_outdir) or args.force.lower() == 'y':
        # Manage all possible genomic variant input files
        if args.vcf is not None:
            report.write(time_format + 'Reading and splitting input file. \n')
            logger.info('Reading and splitting input file.')
            # for loop in case we have multiple inputs to read from a list of files
            for f in args.vcf:
                # check if input is a file
                if isfile(f) == 'list_files':
                    with open(f) as list_var_files:
                        var_f = list_var_files.read().splitlines()
                        logger.info(
                            'Input variants file contains a list of variants files to process.')
                        # for every prot id
                        for var_infile in var_f:
                            # detect the format of the vcf file(s), either .vcf or .vep
                            input_format = detect_format(var_infile)
                            # If vcf transform into vep format and split
                            if input_format == "vcf":
                                # change input format if file doesn't exists or overwrite is True
                                if os.path.isfile(out_file) is False or args.force.lower() == 'y':
                                    # log info
                                    logger.info(
                                        'Input file' + var_infile + ' is in .vcf format.')
                                    # split vcf file
                                    varfile.vcf(var_infile, out_dir,
                                                out_file, vardb_outdir, args.force, log_dir, args.parallel)
                                    # log info
                                    logger.info(
                                        var_infile + ' has been splitted.')
                            # If vep, only split
                            elif input_format == "vep":
                                # split if empty dir or overwrite is True
                                if not os.listdir(vardb_outdir) or args.force.lower() == 'y':
                                     # log info
                                    logger.info(
                                        'Input file' + var_infile + ' is in .vep format.')
                                    # split vep file by protein id to speed up the
                                    # mapping process
                                    varfile.vep(
                                        var_infile, vardb_outdir, args.force, log_dir, args.parallel)
                                    # log info
                                    logger.info(
                                        var_infile + ' has been splitted.')
                            else:
                                # log info
                                spinner.warn('Warning: input file', var_infile,
                                             'is neither in vep nor vcf format.')
                                logger.warning(' Input file', var_infile,
                                               'is not in vep nor vcf format.')
                                continue

                        report.write(
                            time_format + 'Splitting process done.\n')

                elif isfile(f) == 'is_file':
                    # change input format if file doesn't exists or overwrite is True
                    if not os.listdir(vardb_outdir) or args.force.lower() == 'y':
                        logger.info(
                            'Input file is a single file. Starting splitting process.')

                        # detect the format of the vcf file(s), either .vcf or .vep
                        input_format = detect_format(f)

                        # If vcf transform into vep format and split
                        if input_format == "vcf":
                            # change input format if file doesn't exists or overwrite is True
                            if os.path.isfile(out_file) is False or args.force.lower() == 'y':
                                # log info
                                logger.info('Input file' + f +
                                            ' is in .vcf format.')
                                # split vcf file
                                varfile.vcf(f, out_dir,
                                            out_file, vardb_outdir, args.force, log_dir, args.parallel)
                                # log info
                                logger.info(
                                    f + ' has been splitted.')
                                report.write(
                                    time_format + 'Splitting process done.\n')

                        # If vep, only split
                        elif input_format == "vep":
                            print("HOLA K ASE")
                            # split if empty dir or overwrite is True
                            if not os.listdir(vardb_outdir) or args.force.lower() == 'y':
                                # log info
                                logger.info('Input file' + f +
                                            ' is in .vcf format.')
                                # split vep file by protein id to speed up the
                                # mapping process
                                varfile.vep(f, vardb_outdir,
                                            args.force, log_dir, args.parallel)
                                # log info
                                logger.info(
                                    f + ' has been splitted.')
                                report.write(
                                    time_format + 'Splitting process done.\n')

                        else:
                            # log info
                            spinner.warn('Warning: input file', f,
                                         'is neither in vep nor vcf format.')
                            logger.warning('Input file', f,
                                           'is not in vep nor vcf format.')

                            continue
                else:
                    spinner.warn(
                        'The input is neither a file(s) or a file containing a list of files')
                    logger.error(
                        'The input is neither a file(s) or a file containing a list of files')

        # If MAF transform into VEP format and split
        elif args.maf is not None:
            report.write(time_format + 'Reading and splitting input file. \n')
            logger.info('Reading and splitting input file.')
            # for loop in case we have multiple inputs to read from a list of files
            for f in args.maf:
             # check if input is a file
                if isfile(f) == 'list_files':
                    # try:
                    with open(f) as list_var_files:
                        var_f = list_var_files.read().splitlines()
                        logger.info(
                            'Input variants file contains a list of variants files to process.')
                        # for every prot id
                        for var_infile in var_f:
                            # change input format if file doesn't exists or overwrite is True
                            if os.path.isfile(out_file) is False or args.force.lower() == 'y':
                                # split MAF file
                                varfile.maf(var_infile, out_dir,
                                            out_file, vardb_outdir, args.force, log_dir, args.parallel)
                                # log info
                                logger.info('File has been splitted.')
                                report.write(
                                    time_format + 'Splitting process done.\n')

                            else:
                                logger.warning(
                                    vardb_outdir + ' already exists and force overwrite option is not active.')

                elif isfile(f) == 'is_file':
                    # change input format if file doesn't exists or overwrite is True
                    if not os.listdir(vardb_outdir) or args.force.lower() == 'y':
                        # split MAF file
                        varfile.maf(f, out_dir,
                                    out_file, vardb_outdir, args.force, log_dir, args.parallel)
                        # log info
                        logger.info('File has been splitted.')
                        report.write(
                            time_format + 'Input file is in .maf format. Splitting process done.\n')
                    else:
                        logger.warning(
                            vardb_outdir + ' already exists and force overwrite option is not active.')

                else:
                    spinner.warn(
                        'The input is neither a file(s) or a file containing a list of files')
                    logger.error(
                        'The input is neither a file(s) or a file containing a list of files')

        elif args.vep is not None:
            spinner.warn(text='VEP option will be available soon. Using VarMap db instead. \
    Otherwise, please provide your own vcf file with the -vcf option.\n')
            try:
                run_vep()
            except:
                vardb_outdir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'
                spinner.info('Using VarMap db\n')
                exit(-1)
        elif args.varmap is not None:
            spinner.info('Using VarMap db')
            vardb_outdir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'
    else:
        print('A variants database already exists.')

    # record total time of execution in report file
    end = time.time()
    report.write(
        time_format + 'Generation of genomic variants DB in ' +
        vardb_outdir + ' took ' + str(round(end-start, 2)) + 's\n')
    report.close()
    # print in console result
    spinner.stop_and_persist(symbol='\U0001F4CD',
                             text=' makevariantsdb is done. in ' +
                             str(round(end-start, 2)) + 's.')
