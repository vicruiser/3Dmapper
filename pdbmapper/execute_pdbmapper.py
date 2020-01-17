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

    # set up the logging
    logger = open(os.path.join(args.out, 'pdbmapper.log'), 'w')
    logger.write(description)
    logger.write(epilog)
    logger.write('''
    Command line input: 
    -------------------
    \n''')
    logger.write((" ".join(sys.argv)) + '\n' + '\n' + '\n')
    time_format = '[' + time.ctime(time.time()) + '] '
    # create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        spinner.info(text="Directory " + args.out + " created.\n")
        logger.write(time_format + "Directory " + args.out + ' created. \n')
    else:
        spinner.info(
            text=args.out + " is an existing directory. Results will be written in there.\n")
        logger.write(time_format + args.out +
                     'is an existing directory. Results will be written in there. \n')
    # create chimera scripts:
    if args.chimera is not None:
        # chimera()
        pass

    # run PDBmapper
    if args.protid:
        # create output dir if it doesn't exist
        os.makedirs(args.out, exist_ok=True)
        # compute total time of running PDBmapper
        start = timer()
        # decorator to monitor function

        @tags(text_start="Running PDBmapper...",
              text_succeed=" Running PDBmapper...done.",
              text_fail=" Running PDBmapper...failed!",
              emoji=DNA)
        # define function to run PDBmapper with the decorator
        def f():
            # PDBmapper accepts single or multiple protein ids
            # as input as well as prot ids stored in a file
            for protids in args.protid:
                # check if input is a file
                try:
                    with open(protids) as f:
                        lines = f.read().splitlines()
                        # set variable input as file
                        input = "file"

                except:
                    # set variable input as not file
                    input = "not_file"

                if input == "file":
                    for protid in lines:
                        try:
                            # for pids in lines:
                            ensemblIDs = translate_ensembl(
                                protid, args.filter_iso)
                            geneid = ensemblIDs['geneID']
                            transcriptID = ensemblIDs['transcriptID']
                        except IOError:
                            log = open(os.path.join(
                                args.out, 'log_ensembl.File'), 'a')
                            log.write('Warning: ' + protid +
                                      ' has no ENGS.\n')
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
                            continue

                # input is not a file but one or more protein ids
                # given in command line
                elif input == "not_file":
                    # for prot id get the gene id
                    protid = protids
                    try:
                        ensemblIDs = translate_ensembl(
                            protid, args.filter_iso)
                        geneid = ensemblIDs['geneID']
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
                            next
                    except IOError:
                        next
                else:
                    print("wrong input!!")
        # execute function
        f()
        # time execution
        end = timer()
        finish = end - start
        finish = round(finish/60, 3)
        # print result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text='Congratulations!. PDBmapper has run in ' +
                                 str(finish) + ' minutes.')

        log_file.write(sttime + ('Congratulations!. PDBmapper has run in ',
                                 str(finish), ' minutes.'))
        log_file.close()
