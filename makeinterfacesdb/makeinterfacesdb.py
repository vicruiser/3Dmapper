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


class generateIntDB:

    def stats(int_infile, intdb_outdir):
        # count after transforming to vep format the total number
        # of input interfaces, number of genes, etc

        # number of interfaces
        n_int_cmd = "awk '{{print $1,$2,$3,$4,$10}}' {} | uniq | wc -l"
        n_int, err1 = call_subprocess(
            n_int_cmd).format(int_infile)

        # number of proteins: count total number of splitted files
        n_prot_cmd = "wc -l {}"
        n_prot, err2 = call_subprocess(n_prot_cmd(intdb_outdir))

        # number of interfaces ligand
        n_int_ligand_cmd = "awk '{{print $1,$2,$3,$4,$10}}' | grep 'ligand' | uniq | wc -l"
        n_int_ligand, err3 = call_subprocess(
            n_int_ligand_cmd).format(int_infile)

        # number of interfaces protein
        n_int_prot_cmd = "awk '{{print $1,$2,$3,$4,$10}}' | grep 'protein' | uniq | wc -l"
        n_int_prot, err3 = call_subprocess(
            n_int_prot_cmd).format(int_infile)

        # number of interfaces nucleic
        n_int_nucleic_cmd = "awk '{{print $1,$2,$3,$4,$10}}' | grep 'nucleic' | uniq | wc -l"
        n_int_nucleic, err3 = call_subprocess(
            n_int_nucleic_cmd).format(int_infile)

        n_int, n_prot, n_int_ligand, n_int_prot, n_int_nucleic =
        n_int.decode('utf-8'), n_prot.decode('utf-8'),
        n_int_ligand.decode('utf-8'), n_int_prot.decode('utf-8'),
        n_int_nucleic.decode('utf-8')

        with open(os.path.join(intdb_outdir, 'makevariantsdb_stats.info'), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([n_int, "n_interfaces"])
            writer.writerow([n_prot, "n_ENSP"])
            writer.writerow([n_int_ligad, "n_interfaces_ligand"])
            writer.writerow([n_int_prot, "n_interfaces_protein"])
            writer.writerow([n_int_nucleic, "n_interfaces_nucleic"])

        stats_message = ('''

        Stats
        -----
         - Total number of input interfaces ids: {}
         - Total number of corresponding proteins (total number of splitted files): {}
         - Total number of interfaces of type ligand: {}
         - Total number of interfaces of type protein: {}
         - Total number of interfaces of type nucleic: {}
        ''').format(str(n_int), str(n_prot), str(n_int_ligand), str(n_int_prot), str(n_int_nucleic))

        return stats_message


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
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")
    # parse command line options
    args = parse_commandline()

    # initialize class
    makedb = generateIntDB()

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

        if args.force is True:
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
                            stats_message = intdb_file.stats(
                                int_infile, intdb_outdir)
                            # log info
                            logger.info(
                                int_infile + ' has been splitted successfully.')

                elif isfile(f) == 'is_file':
                    # split interface db
                    split('ENSP', f, intdb_outdir,
                          'txt', args.force, log_dir)
                    stats_message = intdb_file.stats(
                        f, intdb_outdir)
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
                    time_format + 'Generation of interfaces DB in ' + intdb_outdir + ' took ' +
                    str(datetime.timedelta(seconds=end-start)) + 's\n')
                report.write(stats_message)
                report.close()
        else:
            makedb.log(
                'A variants database already exists. Not overwritting files.')
            spinner.stop_and_persist(symbol='\U0001F4CD',
                                     text=' A variants database already exists. Not overwritting files.')
            report.close()
