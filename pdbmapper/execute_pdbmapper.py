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
import multiprocessing as mp

from halo import Halo
from timeit import default_timer as timer
from subprocess import call

# import functions from scripts
from .parse_argv import parse_commandline
from .translate_ensembl import translate_ensembl
from .PDBmapper import PDBmapper
from .decorator import tags
from .logger import get_logger
from .input_isfile import isfile
from .run_subprocess import call_subprocess
from .pdbmapper_wrapper import wrapper

pd.options.mode.chained_assignment = None
#print("Number of processors: ", mp.cpu_count())

# define main function


class pdbmapper:

    def run():
    def run_in_parallel():

    def stats():
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

    '''

    epilog = \
        '''
          -------------------------------------------------------------------------
         |  Copyright (c) 2019 Victoria Ruiz --                                    |
         |  vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          -------------------------------------------------------------------------

        '''

    # Emojis
    DNA = '\U0001F9EC'
    # spinner
    spinner = Halo(text='Loading', spinner='dots12', color="red")
    # parse command line options
    args = parse_commandline()

    # print ascii art
    print(description)
    print(epilog)
    # time message
    time_format = '[' + time.ctime(time.time()) + '] '

    # create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        spinner.info(text="Directory " + args.out + " created.\n")
        out_message = "Directory " + args.out + ' created.'

    else:
        spinner.info(
            text=args.out + " is an existing directory. Results will be written in there.\n")
        out_message = args.out + ' is an existing directory. Results will be written in there.'

  # set up the logging
    logger = get_logger('main', args.out)
    logger.info(out_message)

    # set up the results report
    report = open(os.path.join(args.out, 'pdbmapper.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('Command line input:\n' +
                 '-------------------\n' +
                 '\n')
    report.write((" ".join(sys.argv)) + '\n' + '\n' + '\n')
    report.write(time_format + out_message + '\n')

    # if parallel is True
    if args.parallel is True:
        pool = mp.Pool(mp.cpu_count()-1)

        # results = pool.map(howmany_within_range_rowonly, [row for row in data])
        # pool.close()

   # # create chimera scripts:
   # if args.chimera is not None:
   #     # chimera()
   #     pass

   # run PDBmapper
    if args.varid:
        # execute main function and compute executiong time
        start = time.time()
        logger.info('Running PDBmapper...')
        report.write(time_format + 'Running PDBmapper...\n')
        spinner.start(text=' Running PDBmapper...')
        # find variants index file
        index_file = glob.glob(os.path.join(args.vardb, '*.index'))[0]

        for ids in args.varid:
            # run PDBmapper
            if isfile(ids) == 'yes':
                with open(ids) as list_varids:
                    for id in list_varids:
                        # remove \n from the end
                        id = id.replace('\n', '')
                        # grep variant in index file created with makevariantsdb
                        cmd = ('grep \'{}\' {}').format(id, index_file)
                        # call subprocess
                        out, err = call_subprocess(cmd)
                        if err is None and out != b'':
                            toprocess = out.decode('utf-8')
                            # geneid correspond to the 2nd column in the line
                            transcriptid = toprocess.split(" ")[2].strip()
                            # execute PDBmapper
                            wrapper(transcriptid,
                                    args.intdb,
                                    args.vardb,
                                    args.out,
                                    args.pident,
                                    args.consequence,
                                    id)
                        else:
                            logger.error(
                                'Wrong input: {} is not a recognizable variant id'.format(id))
                            continue

            elif isfile(ids) == 'no':
                id = ids
                # grep variant in index file created with makevariantsdb
                cmd = ('grep \'{}\' {}').format(id, index_file)
                # call subprocess
                out, err = call_subprocess(cmd)
                if err is None and out != b'':
                    toprocess = out.decode('utf-8')
                    # geneid correspond to the 2nd column in the line
                    transcriptid = toprocess.split(" ")[2].strip()
                    # execute PDBmapper
                    wrapper(transcriptid,
                            args.intdb,
                            args.vardb,
                            args.out,
                            args.pident,
                            args.consequence,
                            id)
                else:
                    logger.error(
                        'Wrong input: {} is not a recognizable variant id'.format(id))
                    continue
            else:
                logger.error(
                    'The input variants ids provided are not in a valid format.')
                spinner.fail(" Running PDBmapper...failed!")
                report.write(time_format + " Running PDBmapper...failed!")
                raise IOError

        # execute main function and compute executiong time
        end = time.time()
        logger.info('Done.')
        report.write(time_format + 'Congratulations!. PDBmapper has run in ' +
                     str(round(end-start, 2)) + 's.')
        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text='Congratulations!. PDBmapper has run in ' +
                                 str(round(end-start, 2)) + 's.')

    if args.ensemblid:
        # execute main function and compute executiong time
        start = time.time()
        logger.info('Running PDBmapper...')
        report.write(time_format + 'Running PDBmapper...\n')
        spinner.start(text=' Running PDBmapper...')
        # PDBmapper accepts single or multiple protein ids
        # as input as well as prot ids stored in a file
        for ids in args.ensemblid:
            # check if input is a file
            if isfile(ids) == "yes":
                with open(ids) as list_ensemblids:
                    logger.info(
                        'Input variants file contains a list of ensembl ids to process.')
                    # for every ensembl id
                    for ensemblid in list_ensemblids:
                        # remove \n from the end
                        ensemblid = ensemblid.replace('\n', '')
                        # run PDBmapper
                        try:
                            wrapper(ensemblid,
                                    args.intdb,
                                    args.vardb,
                                    args.out,
                                    args.pident,
                                    args.consequence)
                        except:
                            continue

            # input is not a file but one or more protein ids
            # given in command line
            elif isfile(ids) == "no":
                # for prot id get the gene id
                ensemblid = ids
                # run PDBmapper
                wrapper(ensemblid,
                        args.intdb,
                        args.vardb,
                        args.out,
                        args.pident,
                        args.consequence)
                continue
            else:
                logger.error(
                    'The input Ensembl ids provided are not in a valid format.')
                spinner.fail(" Running PDBmapper...failed!")
                report.write(time_format + " Running PDBmapper...failed!")
                raise IOError

        # compute statistics
        # pdbmapper.stats()

        # Compute execution time
        end = time.time()
        logger.info('Done.')
        report.write(time_format + 'Congratulations!. PDBmapper has run in ' +
                     str(round(end-start, 2)) + 's.')
        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text='Congratulations!. PDBmapper has run in ' +
                                 str(round(end-start, 2)) + 's.')
