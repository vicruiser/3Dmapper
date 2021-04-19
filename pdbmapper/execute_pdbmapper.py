# coding: utf-8

# Import necesary modules
from .stats import stats
from .pdbmapper_wrapper import wrapper
from .run_subprocess import call_subprocess
from .input_isfile import isfile
from .logger import get_logger
from .decorator import tags
from .PDBmapper import PDBmapper
from .translate_ensembl import translate_ensembl
from .parse_argv import parse_commandline
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
import shutil


import pandas as pd
import numpy as np
import multiprocessing as mp

from halo import Halo
from timeit import default_timer as timer
from subprocess import call
from tabulate import tabulate
from joblib import Parallel, delayed, parallel_backend


pd.options.mode.chained_assignment = None

# define main function


class MapTools:

    time = '[' + time.ctime(time.time()) + '] '

    def log(self, message, report, logger):

        logger.info(message)
        report.write(self.time + message + '\n')

    def run():
        pass

    def run_in_parallel():
        pass


def finish_message(logger, report, time_format, start, spinner):
    end = time.time()
    logger.info('Done.')
    report.write(time_format + 'PDBmapper has finished successfully. Total time: ' +
                 str(datetime.timedelta(
                     seconds=round(end-start))) + '\n \n')
    report.write('Stats\n')
    report.write('-----\n')
    # report.write(stats_message)
    report.write('\n')
    # report.write(stats_table)
    report.write('\n')
    # print in console result
    spinner.stop_and_persist(symbol='\U0001F9EC',  # '\U0001F9EC'
                             text=' PDBmapper has finished successfully. Total time: ' +
                             str(datetime.timedelta(
                                 seconds=round(end-start))))


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

    epilog = '''
          ------------------------------------------------------------------------------------
         |  Copyright (c) 2019 Victoria Ruiz --                                               |
         |  victoria.ruizserrra@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          ------------------------------------------------------------------------------------

        '''

    # Emojis
    DNA = '\U0001F9EC'
    # spinner
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")
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

    # if overwrite option is True
    if args.force is True:
        # if os.path.exists(os.path.join(args.out, 'pdbmapper.log')):
        #    os.remove(os.path.join(args.out, 'pdbmapper.log'))
        fileList = glob.glob(os.path.join(args.out, 'setID*.txt'))
        # Iterate over the list of filepaths & remove each file.
        for filePath in fileList:
            try:
                os.remove(filePath)
            except:
                pass
        for item in os.listdir(args.out):
            if item in ['hdf5', 'csv']:
                shutil.rmtree(os.path.join(args.out, item))

    elif args.append is True:
        spinner.info(
                text=' Directory ' + args.out + ' is not empty. Appending files. ')
        
    else:
        for item in os.listdir(args.out):
            if item in ['hdf5', 'csv']:
                logger.warning(
                    'Directory ' + args.out + ' is not empty. Not overwritting files. ' +
                    'Please select option --force or --apend or specify a different output dir.')
                spinner.warn(
                    text=' Directory ' + args.out + ' is not empty. Not overwritting files. ' +
                    'Please select option --force or --apend specify a different output dir.')
                exit(-1)

    out_hdf = args.out + '/hdf5'
    if args.hdf is True:
        if not os.path.exists(out_hdf):
            os.mkdir(out_hdf)
            spinner.info(text="Directory " + out_hdf + " created.\n")
            out_message = "Directory " + out_hdf + ' created.'
        else:
            spinner.info(
                text=out_hdf + " is an existing directory. Results will be written in there.\n")
            out_message = out_hdf + ' is an existing directory. Results will be written in there.'

    out_csv = args.out + '/csv'
    if args.csv is True:
        if not os.path.exists(out_csv):
            os.mkdir(out_csv)
            spinner.info(text="Directory " + out_csv + " created.\n")
            out_message = "Directory " + out_csv + ' created.'
        else:
            spinner.info(
                text=out_csv + " is an existing directory. Results will be written in there.\n")
            out_message = out_csv + ' is an existing directory. Results will be written in there.'

    # set up the results report
    report = open(os.path.join(args.out, 'pdbmapper.report'), 'w')
    report.write(description)
    report.write(epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    report.write(time_format + out_message + '\n')

    # initialize class
    maptools = MapTools()

    # print("Number of processors: ", mp.cpu_count())

    # if parallel is True
    if args.parallel is True:
        if args.njobs:
            num_cores = int(args.njobs)
        else:
            num_cores = mp.cpu_count()-1
    else:
        num_cores = 1

   # # create chimera scripts:
   # if args.chimera is not None:
   #     # chimera()
   #     pass

   # run PDBmapper
    index_file = glob.glob(os.path.join(args.vardb, '*.index'))[0]
    if args.varid:
        # execute main function and compute executiong time
        start = time.time()
        logger.info('Running PDBmapper...')
        report.write(time_format + 'Running PDBmapper...\n')
        if args.verbose:
            spinner.start(text=' Running PDBmapper...')
        else:
            print('Running PDBmapper...')
        # find variants index file

        for ids in args.varid:
           # print(ids)
            # run PDBmapper
            if isfile(ids) == 'yes':
                input = 'file'
                with open(ids) as list_varids:
                    for id in list_varids:
                        # remove \n from the end
                        id = id.replace('\n', '')
                        # grep variant in index file created with makevariantsdb
                        cmd = ('grep \'\\b{}\\b\' {}').format(id, index_file)
                        # call subprocess
                        out, err = call_subprocess(cmd)
                        if err is None and out != b'':
                            toprocess = out.decode('utf-8')
                            # geneid correspond to the 2nd column in the line
                            transcriptid = [s for s in toprocess.split(
                                " ") if 'ENST' in s]
                            #transcriptid = toprocess.split(" ")[2].strip()
                            # execute PDBmapper
                            Parallel(n_jobs=num_cores)(delayed(wrapper)(t,
                                                                        args.psdb,
                                                                        args.vardb,
                                                                        args.out,
                                                                        args.pident,
                                                                        args.evalue,
                                                                        args.isoform,
                                                                        args.consequence,
                                                                        args.loc,
                                                                        index_file,
                                                                        args.uniprot,
                                                                        id,
                                                                        args.csv,
                                                                        args.hdf)
                                                       for t in transcriptid)
                        else:
                            logger.error(
                                'Wrong input: {} is not a recognizable variant id'.format(id))
                            continue
            # given in command line
            elif isfile(ids) == "no":
                #input = 'not_file'
                # break
                cmd = ('grep \'\\b{}\\b\' {}').format(ids, index_file)
                # call subprocess
                out, err = call_subprocess(cmd)
                if err is None and out != b'':
                    toprocess = out.decode('utf-8').split(" ")
                    # geneid correspond to the 2nd column in the line
                    transcriptid = [
                        s for s in toprocess if ('ENST' or '-') in s]
                    if any(transcriptid) is False and len(toprocess) > 0:
                        transcriptid = ['-']
                    Parallel(n_jobs=num_cores)(delayed(wrapper)(t,
                                                                args.psdb,
                                                                args.vardb,
                                                                args.out,
                                                                args.pident,
                                                                args.evalue,
                                                                args.isoform,
                                                                args.consequence,
                                                                args.loc,
                                                                index_file,
                                                                args.uniprot,
                                                                ids,
                                                                args.csv,
                                                                args.hdf)
                                               for t in transcriptid)
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

            finish_message(logger, report, time_format, start, spinner)

            # if input == 'not_file':
            #     cmd = ('grep -E \'{}\' {}').format('|'.join(args.varid), index_file)
            #     # call subprocess
            #     out, err = call_subprocess(cmd)
            #     #print(out, err)
            #     if err is None and out != b'':
            #         toprocess = out.decode('utf-8')
            #         # geneid correspond to the 2nd column in the line
            #         transcriptid = [s for s in toprocess.split(
            #                         " ") if 'ENST' in s]
            #     #print(transcriptid)
            #     Parallel(n_jobs=num_cores)(delayed(wrapper)(transcriptid[i],
            #                                                 args.psdb,
            #                                                 args.vardb,
            #                                                 args.out,
            #                                                 args.pident,
            #                                                 args.isoform,
            #                                                 args.consequence,
            #                                                 args.loc,
            #                                                 args.uniprot,
            #                                                 varids[i],
            #                                                 args.csv,
            #                                                 args.hdf)
            #                             for i in len(transcriptid))

            # if not any(fname.endswith('.File') for fname in os.listdir(args.out)):
            #     logger.warning(
            #         'Error: Input ensembl ids has no mapping variants.')
            #     spinner.warn(
            #         text=' Input ensembl ids has no mapping variants.')
            #     exit(-1)

        # elif isfile(ids) == 'no':
        #         id = ids
        #         # grep variant in index file created with makevariantsdb
        #         cmd = ('grep \'{}\' {}').format(id, index_file)
        #         # call subprocess
        #         out, err = call_subprocess(cmd)
        #         if err is None and out != b'':
        #             toprocess = out.decode('utf-8')
        #             # geneid correspond to the 2nd column in the line
        #             transcriptid = toprocess.split(" ")[2].strip()
        #             # execute PDBmapper
        #             wrapper(transcriptid,
        #                     args.psdb,
        #                     args.vardb,
        #                     args.out,
        #                     args.pident,
        #                     args.consequence,
        #                     id)
        #         else:
        #             logger.error(
        #                 'Wrong input: {} is not a recognizable variant id'.format(id))
        #             continue

        # # execute main function and compute executiong time
        # end = time.time()
        # logger.info('Done.')
        # report.write(time_format + 'Congratulations!. PDBmapper has run in ' +
        #              str(round(end-start, 2)) + 's.')
        # # print in console result
        # spinner.stop_and_persist(symbol='\U0001F4CD',
        #                          text='Congratulations!. PDBmapper has run in ' +
        #                          str(round(end-start, 2)) + 's.')

    if args.protid:
        # execute main function and compute executiong time
        start = time.time()
        logger.info('Running PDBmapper...')
        report.write(time_format + 'Running PDBmapper...\n')
        if args.verbose:
            spinner.start(text=' Running PDBmapper...')
        else:
            print('Running PDBmapper...')
        # PDBmapper accepts single or multiple protein ids
        # as input as well as prot ids stored in a file
        for ids in args.protid:
            # check if input is a file
            if isfile(ids) == "yes":
                input = 'file'
                with open(ids) as list_protids:
                    logger.info(
                        'Input variants file contains a list of ensembl ids to process.')
                    # for every ensembl id
                    Parallel(n_jobs=num_cores)(delayed(wrapper)(protid.replace('\n', ''),
                                                                args.psdb,
                                                                args.vardb,
                                                                args.out,
                                                                args.pident,
                                                                args.evalue,
                                                                args.isoform,
                                                                args.consequence,
                                                                args.loc,
                                                                index_file,
                                                                args.uniprot,
                                                                None,
                                                                args.csv,
                                                                args.hdf)
                                               for protid in list_protids)

            # given in command line
            elif isfile(ids) == "no":
                input = 'not_file'
                break
            elif isfile(ids) == "not_recognized":
                maptools.log('The input is neither an id(s) or a file containing a list of ids.',
                             report, logger)
                spinner.fail(
                    'The input is neither an id(s) or a file containing a list of ids.')
                exit(-1)

        # for prot id get the gene id
        if input == 'not_file':
            # print(args.protid)
            Parallel(n_jobs=num_cores)(delayed(wrapper)(ids,
                                                        args.psdb,
                                                        args.vardb,
                                                        args.out,
                                                        args.pident,
                                                        args.evalue,
                                                        args.isoform,
                                                        args.consequence,
                                                        args.loc,
                                                        index_file,
                                                        args.uniprot,
                                                        None,
                                                        args.csv,
                                                        args.hdf)
                                       for ids in args.protid)

        if not any(fname.endswith('.txt') for fname in os.listdir(args.out)):
            logger.warning(
                'Error: Input ensembl ids has no mapping variants.')
            spinner.warn(
                text=' Input ensembl ids has no mapping variants.')
            exit(-1)

        # var_statsfile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
        #    args.vardb, 'makevariantsdb_stats.info'))[0]))

        # int_statsfile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
        #    args.psdb, 'makeinterfacesdb_stats.info'))[0]))

        # mapped_infofile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
        #    args.out, 'MappedVariants*' + str(args.pident) + '_isoform_' +
        #     '_'.join(args.isoform) + '_consequence_' + '_'.join(args.consequence) + '*File'))[0]))

        # stats_message, stats_table = stats(
        #    var_statsfile, int_statsfile, mapped_infofile, args.out)

        # Compute execution time
        finish_message(logger, report, time_format, start, spinner)
