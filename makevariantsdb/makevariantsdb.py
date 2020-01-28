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
import datetime
import csv

import os.path as path
import pandas as pd
import numpy as np

from halo import Halo
from timeit import default_timer as timer
from subprocess import call
from multiprocessing import Pool
from memory_profiler import memory_usage

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
from .run_subprocess import call_subprocess

# num_cpus = psutil.cpu_count(logical=False)


class generateVarDB:

    time = '[' + time.ctime(time.time()) + '] '

    def log(self, message, report, logger):

        logger.info(message)
        report.write(self.time + message + '\n')

    def stats(self, var_infile, vardb_outdir):
        # count after transforming to vep format the total number
        # of input variants, number of genes, etc

        # number of variants: assuming that they are stored in the first column
        n_variants_cmd = "grep -v 'Gene' {} | awk '{{print $1}}' | uniq | wc -l"
        n_variants, err1 = call_subprocess(n_variants_cmd.format(var_infile))

        # number of genes: count total number of splitted files
        n_genes_cmd = "ls {} | wc -l"
        n_genes, err2 = call_subprocess(n_genes_cmd.format(vardb_outdir))

        # number of trasncripts: assuming that they are stored in the 5th column
        n_features_cmd = "grep -v 'Gene' {} | awk '{{print $5}}' | uniq | wc -l"
        n_features, err3 = call_subprocess(n_features_cmd.format(var_infile))

        n_variants, n_genes, n_features = n_variants.decode(
            'utf-8').rstrip(), n_genes.decode('utf-8').rstrip(), n_features.decode('utf-8').rstrip()

        # store results
        d = {'n_variants': [n_variants],
             'n_genes': [n_genes],
             'n_features': [n_features]}
        df = pd.DataFrame(data=d)
        df.to_csv(os.path.join(vardb_outdir, 'makevariantsdb_stats.info'),
                  sep='\t', encoding='utf-8', index=False)
        stats_message = ('''

        Stats
        -----
         - Total number of input variants: {}
         - Total number of corresponding genes (total number of splitted files): {}
         - Total number of corresponding features: {}
                 ''').format(str(n_variants), str(n_genes), str(n_features))

        return stats_message

    def vcf(self, var_infile, out_dir, out_file, overwrite, log_dir, parallel=False):

        # from vcf to vep
        vcf2vep(var_infile, out_dir,
                out_file, overwrite, log_dir, parallel)
        # add header to resulting vep file
        add_header(out_file)

    def vep(self, var_infile, vardb_outdir, overwrite, log_dir, parallel=False):

        # split vep file by protein id to speed up the
        # mapping process
        split('ENSG', var_infile, vardb_outdir,
              'vep', overwrite, log_dir, parallel)

    def maf(self, var_infile, out_dir, out_file, vardb_outdir, overwrite, log_dir, report, logger, parallel=False):
        # from vcf to vep
        maf2vep(var_infile, out_dir,
                out_file, overwrite, log_dir)
        # split vep file by protein id to speed up the
        # mapping process
        var_infile = out_file

        self.vep(
            var_infile, vardb_outdir, overwrite, log_dir, parallel)

        # logging
        self.log('Splitting process is done. Calculating stats info...',
                 report, logger)
        # report stats
        stats_message = self.stats(
            var_infile, vardb_outdir)
        # log info
        self.log('Calculating stats info...done.',
                 report, logger)

        return stats_message

    def wrapper(self, input_format, var_infile, out, log_dir,
                report, logger, spinner, overwrite=False, parallel=False):

        # created by default
        out_dir = os.path.join(out, 'DBs')
        out_file = os.path.join(
            out_dir, 'variants.vep')  # created by default
        # set output dir to split vep
        vardb_outdir = os.path.join(out_dir, 'varDB')  # created by default

        # If vcf transform into vep format and split
        if input_format == "vcf" or input_format == "vep":
            # change input format if file doesn't exists or overwrite is True
            if input_format == "vcf":
                # split vcf file
                self.vcf(var_infile, out_dir,
                         out_file, overwrite,
                         log_dir, parallel)
                # logging
                self.log('Input vcf file converted to vep format. Splitting vep file...',
                         report, logger)

                var_infile = out_file

            self.vep(
                var_infile, vardb_outdir, overwrite, log_dir, parallel)
            # logging
            self.log('Splitting process is done. Calculating stats info...',
                     report, logger)
            # report stats
            stats_message = self.stats(
                var_infile, vardb_outdir)
            # log info
            self.log('Calculating stats info...done.',
                     report, logger)

            return stats_message

        else:
            # log info
            spinner.warn('Warning: input file' + f +
                         'is neither in vep nor vcf format.')
            makedb.log('Warning: input file' + f +
                       'is neither in vep nor vcf format.',
                       report, logger)
            return IOError


# main function
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
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")
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
    makedb = generateVarDB()

    # set up the results report
    report = open(os.path.join(out_dir, 'makevariantsdb.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('''
        Command line input:
        -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    time_format = '[' + time.ctime(time.time()) + '] '

    # set up a log file
    logger = get_logger('main', out_dir)
    log_dir = out_dir
    # spinner.start(text=' Running makevariantsdb...')
    start = time.time()

    # change input format if file doesn't exists or overwrite is True
    if not os.listdir(vardb_outdir) or args.force is True:
        # Manage all possible genomic variant input files
        if args.vcf is not None:
            # logging
            makedb.log('Starting process of splitting input file of variants.',
                       report, logger)
            # for loop in case we have multiple inputs to read from a list of files
            for f in args.vcf:
                # check if input is a file
                if isfile(f) == 'list_files':
                    with open(f) as list_var_files:
                        # read lines
                        var_f = list_var_files.read().splitlines()
                        logger.info(
                            'Input is a list of files.')
                        # for every prot id
                        for var_infile in var_f:
                            # detect the format of the vcf file(s), either .vcf or .vep
                            input_format = detect_format(var_infile)
                            logger.info(
                                'Input file is in ' + input_format + ' format.')
                            try:
                                stats_message = makedb.wrapper(
                                    input_format, f, args.out, log_dir, report,
                                    logger, spinner, args.force, args.parallel)
                            except IOError:
                                continue

                elif isfile(f) == 'is_file':
                    logger.info(
                        'Input is a single file.')
                    # detect the format of the vcf file(s), either .vcf or .vep
                    input_format = detect_format(f)
                    logger.info(
                        'Input file is in ' + input_format + ' format.')
                    try:
                        stats_message = makedb.wrapper(
                            input_format, f, args.out, log_dir, report,
                            logger, spinner, args.force, args.parallel)
                    except IOError:
                        continue

                elif isfile(f) == 'file_not_recognized':
                    makedb.log('The input is neither a file(s) or a file containing a list of files.',
                               report, logger)
                    spinner.fail(
                        'The input is neither a file(s) or a file containing a list of files')
                    exit(-1)

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
                            'Input is a list of files.')
                        # for every prot id
                        for var_infile in var_f:
                            # split MAF file
                            try:
                                stats_message = makedb.maf(var_infile, out_dir,
                                                           out_file, vardb_outdir,
                                                           args.force, log_dir, report, logger, args.parallel)
                            except IOError:
                                continue

                elif isfile(f) == 'is_file':
                    logger.info(
                        'Input is a single file.')
                    # split MAF file
                    try:
                        stats_message = makedb.maf(f, out_dir,
                                                   out_file, vardb_outdir,
                                                   args.force, log_dir,  report, logger, args.parallel)
                    except IOError:
                        continue

                else:
                    logger.error(
                        'The input is neither a file(s) or a file containing a list of files')
                    spinner.fail(
                        'The input is neither a file(s) or a file containing a list of files')
                    exit(-1)

        # elif args.vep is not None:
        #     spinner.warn(text='VEP option will be available soon. Using VarMap db instead. \
        # Otherwise, please provide your own vcf file with the -vcf option.\n')
        #     try:
        #         run_vep()
        #     except:
        #         vardb_outdir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'
        #         spinner.info('Using VarMap db\n')
        #         exit(-1)
        # elif args.varmap is not None:
        #     spinner.info('Using VarMap db')
        #     vardb_outdir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'

            # record total time of execution in report file
        end = time.time()
        report.write(
            time_format + 'Generation of genomic variants DB in ' +
            vardb_outdir + ' completed successfully. Total time: ' + str(datetime.timedelta(seconds=round(end-start))) + '\n')
        report.write(stats_message)

        # mem_usage = memory_usage(f)
        # print('Memory usage (in chunks of .1 seconds): %s' % mem_usage)
        # print('Maximum memory usage: %s' % max(mem_usage))
        report.close()
        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text=' makevariantsdb process finished. Total time:  ' +
                                 str(datetime.timedelta(
                                     seconds=round(end-start))))
    else:
        makedb.log('A variants database already exists. Not overwritting files.')
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text=' A variants database already exists. Not overwritting files.')
        report.close()
