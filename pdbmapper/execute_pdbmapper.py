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

    # print ascii art
    print(description)
    print(epilog)
    # Emojis
    DNA = '\U0001F9EC'

    # spinner
    spinner = Halo(text='Loading', spinner='dots12', color="red")

    # parse command line options
    args = parse_commandline()

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
    # create chimera scripts:
    if args.chimera is not None:
        # chimera()
        pass

    # run PDBmapper
    if args.varid:
        print(args.vardb)
        index_file = glob.glob(os.path.join(args.vardb, '*.index'))[0]
        print(index_file)
        for ids in args.varid:
            # run PDBmapper
            if isfile(ids) == 'yes':
                for id in ids:

                    cmd = 'grep \'' + id + '\' ' + index_file
                    out, err = call_subprocess(cmd)
                    toprocess = out.decode('utf-8')
            elif isfile(ids) == 'no':
                print(ids)
                cmd = 'grep \'' + ids + '\' ' + index_file
                out, err = call_subprocess(cmd)
                toprocess = out.decode('utf-8')
                geneID = toprocess.split(" ")[1].strip()
                print(j)
            else:
                print("wrong input")

    if args.ensemblid:

        # decorator to monitor function
        @tags(text_start="Running PDBmapper...",
              text_succeed=" Running PDBmapper...done.",
              text_fail=" Running PDBmapper...failed!",
              emoji=DNA)
        # define function to run PDBmapper with the decorator
        def f():
            # PDBmapper accepts single or multiple protein ids
            # as input as well as prot ids stored in a file
            for ids in args.ensemblid:
                # check if input is a file
                if isfile(ids) == "yes":
                    for ensemblid in lines:
                        try:
                            # for pids in lines:
                            ensemblIDs = translate_ensembl(
                                ensemblid, args.filter_iso, args.out)
                            geneid = ensemblIDs['geneID']
                            protid = ensemblIDs['protID']
                            transcriptID = ensemblIDs['transcriptID']
                        except IOError:
                            logger.error('Warning: ' + ensemblid +
                                         ' has no corresponding translation.')
                            continue
                        # run PDBmapper
                        try:
                            PDBmapper(protid,
                                      geneid,
                                      transcriptID,
                                      args.intdb,
                                      args.vardb,
                                      args.out,
                                      args.pident,
                                      args.filter_var)
                        # error handling
                        except IOError:
                            logger.error('Warning: ' + ensemblid +
                                         ' has no mapping variants.')
                            continue

                # input is not a file but one or more protein ids
                # given in command line
                elif isfile(ids) == "no":
                    # for prot id get the gene id
                    ensemblid = ids
                    try:
                        ensemblIDs = translate_ensembl(
                            ensemblid, args.filter_iso, args.out)
                        geneid = ensemblIDs['geneID']
                        protid = ensemblIDs['protID']
                        transcriptID = ensemblIDs['transcriptID']
                        # run PDBmapper
                        try:
                            PDBmapper(protid,
                                      geneid,
                                      transcriptID,
                                      args.intdb,
                                      args.vardb,
                                      args.out,
                                      args.pident,
                                      args.filter_var)
                    # error handling
                        except IOError:
                            logger.error('Warning: ' + ensemblid +
                                         ' has no mapping variants.')
                            next
                    except IOError:
                        logger.error('Warning: ' + ensemblid +
                                     ' has no corresponding translation.')
                        next
                else:
                    logger.error('Wrong input!.')

        # execute main function and compute executiong time
        logger.info('Running PDBmapper...')
        report.write(time_format + 'Running PDBmapper...\n')

        start = time.time()
        f()
        end = time.time()

        logger.info('Done.')
        report.write(time_format + 'Congratulations!. PDBmapper has run in ' +
                     str(round(end-start, 2)) + 's.')

        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text='Congratulations!. PDBmapper has run in ' +
                                 str(round(end-start, 2)) + 's.')
