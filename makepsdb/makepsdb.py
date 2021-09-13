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
import datetime
import csv
import shutil

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
from .input_isfile import isfile
from .run_subprocess import call_subprocess


class generateIntDB:

    time = '[' + time.ctime(time.time()) + '] '

    def log(self, message, report, logger):

        logger.info(message)
        report.write(self.time + message + '\n')

    # def stats(self, int_infile, psdb_outdir):
    #     # count after transforming to vep format the total number
    #     # of input interfaces, number of genes, etc

    #     # number of interfaces
    #     n_int_cmd = " awk 'NR>1 {{print $1,$2,$3,$4,$10}}' {} | wc -l"
    #     n_int, err1 = call_subprocess(
    #         n_int_cmd.format(int_infile))

    #     # number of proteins: count total number of splitted files
    #     n_prot_cmd = "find {} -type f -iname 'ENSP*' | wc -l"
    #     n_prot, err2 = call_subprocess(n_prot_cmd.format(psdb_outdir))

    #     # number of interfaces ligand
    #     n_int_ligand_cmd = "awk '{{print $1,$2,$3,$4,$10}}' {} | grep 'ligand' |  wc -l"
    #     n_int_ligand, err3 = call_subprocess(
    #         n_int_ligand_cmd.format(int_infile))

    #     # number of interfaces protein
    #     n_int_prot_cmd = "awk '{{print $1,$2,$3,$4,$10}}' {} | grep 'protein' |  wc -l"
    #     n_int_prot, err3 = call_subprocess(
    #         n_int_prot_cmd.format(int_infile))

    #     # number of interfaces nucleic
    #     n_int_nucleic_cmd = "awk '{{print $1,$2,$3,$4,$10}}' {} | grep 'nucleic' |  wc -l"
    #     n_int_nucleic, err3 = call_subprocess(
    #         n_int_nucleic_cmd.format(int_infile))

    #     n_int, n_prot, n_int_ligand, n_int_prot, n_int_nucleic = n_int.decode('utf-8').rstrip(), n_prot.decode('utf-8').rstrip(
    #     ), n_int_ligand.decode('utf-8').rstrip(), n_int_prot.decode('utf-8').rstrip(), n_int_nucleic.decode('utf-8').rstrip()

    #     d = {'n_interfaces': [n_int],
    #          'n_ENSP': [n_prot],
    #          'n_interfaces_ligand': [n_int_ligand],
    #          'n_interfaces_protein': [n_int_prot],
    #          'n_interfaces_nucleic': [n_int_nucleic]}
    #     df = pd.DataFrame(data=d)
    #     df.to_csv(os.path.join(psdb_outdir, 'makeinterfacesdb_stats.info'),
    #               sep='\t', encoding='utf-8', index=False)

    #     stats_message = ('''

    #     Stats
    #     -----
    #      - Total number of input interfaces ids: {}
    #      - Total number of corresponding proteins (total number of splitted files): {}
    #      - Total number of interfaces of type ligand: {}
    #      - Total number of interfaces of type protein: {}
    #      - Total number of interfaces of type nucleic: {}
    #     ''').format(str(n_int), str(n_prot), str(n_int_ligand), str(n_int_prot), str(n_int_nucleic))

    #     return stats_message


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

          ---------------  Map genomic variants to protein data in 3D. -----------------

    \n'''

    epilog = \
        '''
          ------------------------------------------------------------------------------------
         |  Copyright (c) 2019 Victoria Ruiz --                                               |
         |  victoria.ruizserrra@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          ------------------------------------------------------------------------------------

        '''
    # print ascii art
    print(description)
    print(epilog)
    # initialize spinner decorator
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")
    # parse command line options
    args = parse_commandline()

    # initialize class
    makedb = generateIntDB()

    # set out dir and out file names
    # created by default
    out_dir = os.path.join(args.out, 'DBs')
    # set output dir to split vep
    psdb_outdir = os.path.join(out_dir, 'psdb')  # created by default
    # create output dir if it doesn't exist
    if args.force is True:
        if os.path.exists(psdb_outdir):
            shutil.rmtree(psdb_outdir)
        os.makedirs(psdb_outdir, exist_ok=True)
    else:
        if os.path.exists(psdb_outdir):
            spinner.warn(
                    text=' Directory ' + psdb_outdir + ' is not empty. Not overwritting files. ' +
                    'Please select option --force or specify a different output dir.')
            exit(-1)
        else: 
            # create output dir if it doesn't exist
            os.makedirs(psdb_outdir, exist_ok=True)
    

    # set up a log file
    logger = get_logger('main', out_dir)
    log_dir = out_dir
    # set up the report
    report = open(os.path.join(out_dir, 'makepsdb.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    time_format = '[' + time.ctime(time.time()) + '] '
    start = time.time()

    if args.psdb is not None:

        if not os.listdir(psdb_outdir) or args.force is True:

            logger.info('Reading and splitting input file.')

            # report info
            report.write(time_format + 'Reading and splitting input file. \n')
            # for loop in case we have multiple inputs to read from a list of files
            for f in args.psdb:
                # check if input is a file
                if isfile(f) == 'list_files':
                    with open(f) as list_int_files:
                        int_f = list_int_files.read().splitlines()
                        logger.info(
                            'Input file contains a list of protein files to process.')
                        # for every prot id
                        for int_infile in int_f:
                            # split interface db
                            split('Protein_accession', int_infile, psdb_outdir,
                                  'txt', args.force, log_dir)
                            # stats_message = makedb.stats(
                            #    int_infile, psdb_outdir)
                            # log info
                            logger.info(
                                int_infile + ' has been splitted successfully.')

                elif isfile(f) == 'is_file':
                    # split interface db
                    split('Protein_accession', f, psdb_outdir,
                          'txt', args.force, log_dir)
                    # stats_message = makedb.stats(
                    #    f, psdb_outdir)
                    # log info
                    logger.info(
                        f + ' has been splitted successfully.')
                elif isfile(f) == 'file_not_recognized':
                    logger.error(
                        'Error: Not such file: \'' + f + '\'')
                    print('Error: Not such file: \'' + f + '\'')
                    exit(-1)

                # finish report
                end = time.time()
                report.write(
                    time_format + 'Reading and splitting input file...done. \n')
                report.write(
                    time_format + 'Generation of protein structures DB in ' + psdb_outdir + ' done. Total time: ' +
                    str(datetime.timedelta(seconds=round(end-start))) + '\n')
                # report.write(stats_message)
                report.close()
                # print in console result
                spinner.stop_and_persist(symbol='\U0001F4CD',
                                         text=' makepsdb process finished. Total time: ' +
                                         str(datetime.timedelta(
                                             seconds=round(end-start))))
        else:
            makedb.log(
                'A protein structures DB already exists. Not overwritting files.', report, logger)
            spinner.stop_and_persist(symbol='\U0001F4CD',
                                     text=' A protein structures DB exists. Not overwritting files.')
            report.close()
