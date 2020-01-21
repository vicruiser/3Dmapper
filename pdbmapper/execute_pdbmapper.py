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
from .translate_ensembl import translate_ensembl
from .PDBmapper import PDBmapper
from .decorator import tags
from .logger import get_logger
from .input_isfile import isfile
from .run_subprocess import call_subprocess
from .pdbmapper_wrapper import wrapper

pd.options.mode.chained_assignment = None


# define main function


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
                                    args.feature,
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
                            args.feature,
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
                        wrapper(ensemblid,
                                args.intdb,
                                args.vardb,
                                args.out,
                                args.pident,
                                args.feature)
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
                        args.feature)
                continue
            else:
                logger.error(
                    'The input Ensembl ids provided are not in a valid format.')
                spinner.fail(" Running PDBmapper...failed!")
                report.write(time_format + " Running PDBmapper...failed!")
                raise IOError
        # Compute execution time
        end = time.time()
        logger.info('Done.')
        report.write(time_format + 'Congratulations!. PDBmapper has run in ' +
                     str(round(end-start, 2)) + 's.')
        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text='Congratulations!. PDBmapper has run in ' +
                                 str(round(end-start, 2)) + 's.')
