# coding: utf-8

# Import necesary modules
from .stats import stats
from .mapper_wrapper import wrapper
from .run_subprocess import call_subprocess
from .input_isfile import isfile
from .logger import get_logger
from .decorator import tags
from .mapper import mapper
from .translate import translate
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
    report.write(time_format + '3Dmapper has finished successfully. Execution time: ' +
                 str(datetime.timedelta(
                     seconds=round(end-start))) + '\n \n')
    spinner.stop_and_persist(symbol='\U0001F9EC',  # '\U0001F9EC'
                             text=' 3Dmapper has finished successfully. Execution time: ' +
                             str(datetime.timedelta(
                                 seconds=round(end-start))))

def aesthetis():
    
    description = '''
    ----------------------------------------------------
    ____  _____                                        
   |___ \|  __ \                                       
     __) | |  | |_ __ ___   __ _ _ __  _ __   ___ _ __ 
    |__ <| |  | | '_ ` _ \ / _` | '_ \| '_ \ / _ \ '__|
    ___) | |__| | | | | | | (_| | |_) | |_) |  __/ |   
   |____/|_____/|_| |_| |_|\__,_| .__/| .__/ \___|_|   
                                | |   | |              
                                |_|   |_|              
    ----------------------------------------------------
    '''

    epilog = '''
          -----------------------------------
         |  >>>Publication link<<<<<         |
         |  victoria.ruiz.serra@gmail.com    |
          -----------------------------------
        '''
    # print ascii art
    print(description)
    print(epilog)


def out_file(out_dir, spinner):
    # create output directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        spinner.info(text="Directory " + out_dir + " created.\n")
        out_message = "Directory " + out_dir+ ' created.'
    else:
        spinner.info(
            text=out_dir + " is an existing directory. Results will be written in there.\n")
        out_message = out_dir + ' is an existing directory. Results will be written in there.'
    return(out_message)



def dest_results(spinner, logger):
    args = parse_commandline()
    if args.force is True:
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
                    'Please select option --force or --append or specify a different output dir.')
                spinner.warn(
                    text=' Directory ' + args.out + ' is not empty. Not overwritting files. ' +
                    'Please select option --force or --append specify a different output dir.')
                exit(-1)

def result_format(arg, out, f_format, spinner, logger):

    out_dest = os.path.join(out, f_format)
    if arg is True:
        if not os.path.exists(out_dest):
            os.makedirs(out_dest, exist_ok=True )
            spinner.info(text="Directory " + out_dest  + " created.\n")
            out_message = "Directory " + out_dest + ' created.'
        else:
            spinner.info(
                text=out_dest  + " is an existing directory. Results will be written in there.\n")
            out_message = out_dest  + ' is an existing directory. Results will be written in there.'
        logger.info(out_message)

def report(out, time_format):
    # set up the results report
    report = open(os.path.join(out, '3dmapper.report'), 'a')
    #report.write(description)
    #report.write(epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    report.write(time_format + out_message + '\n')


def parallel(parallel, njobs): 
    if parallel is True:
        if njobs != 1:
            num_cores = int(njobs)
        else:
            num_cores = mp.cpu_count()-1
    else:
        num_cores = njobs
    return(num_cores)

def start_spinner(verbose, logger, time_format, spinner):
    start = time.time()
    logger.info('Running 3Dmapper...')
    if verbose:
        spinner.start(text=' Running 3Dmapper...')
    else:
        print('Running 3Dmapper...')
    return(start)


def main():
    # parse command line options
    args = parse_commandline()

    # Emojis & spinner
    DNA = '\U0001F9EC'
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")

    # set up the logging
    out_message = out_file(args.out, spinner)
    logger = get_logger('main', args.out)
    logger.info(out_message)
    # print ascii art
    aesthetis()
    # time message
    time_format = '[' + time.ctime(time.time()) + '] '

    dest_results(spinner, logger)

    result_format(args.hdf, args.out, 'hdf5', spinner, logger)
    result_format(args.csv, args.out, 'csv', spinner, logger)

    # set up the results report
    report = open(os.path.join(args.out, '3dmapper.report'), 'a')
    #report.write(description)
    #report.write(epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    report.write(time_format + out_message + '\n')

    # initialize class
    maptools = MapTools()

    # if parallel is True
    num_cores  = parallel(args.parallel, args.njobs)
    # find variants index file created with makevariantsdb
    index_file = glob.glob(os.path.join(args.vardb, '*.index'))[0]
    

    start= start_spinner(args.verbose, logger, time_format, spinner)
    if args.varid:
        # find positions index file
        for ids in args.varid:
            # run PDBmapper
            if isfile(ids) == 'yes':
                input = 'file'
                with open(ids) as list_varids:
                    for id in list_varids:
                        # remove \n from the end
                        id = id.replace('\n', '')
                        # grep position in index file created with makepositionsdb
                        cmd = ('grep \'\\b{}\\b\' {}').format(id, index_file)
                        out, err = call_subprocess(cmd)
                        if err is None and out != b'':
                            toprocess = out.decode('utf-8')
                            transcriptid = [toprocess[2].strip()]

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
                                                                        args.dict_geneprot,
                                                                        #args.uniprot,
                                                                        id,
                                                                        args.csv,
                                                                        args.hdf)
                                                       for t in transcriptid)
                        else:
                            logger.error(
                                'Wrong input: {} is not a recognizable position id'.format(id))
                            continue

            elif isfile(ids) == "no":
                cmd = ('grep \'\\b{}\\b\' {}').format(ids, index_file)
                # call subprocess
                out, err = call_subprocess(cmd)
                if err is None and out != b'':
                    toprocess = out.decode('utf-8').split(" ")
                    transcriptid = [toprocess[2].strip()]
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
                                                                args.dict_geneprot,
                                                               # args.uniprot,
                                                                ids,
                                                                args.csv,
                                                                args.hdf)
                                               for t in transcriptid)
                else:
                    logger.error(
                        'Wrong input: {} is not a recognizable position id'.format(ids))
                    continue

            else:
                logger.error(
                    'The input positions ids provided are not in a valid format.')
                spinner.fail(" Running 3Dmapper...failed!")
                report.write(time_format + " Running 3Dmapper...failed!")
                raise IOError
            finish_message(logger, report, time_format, start, spinner)

    if args.prot_id:
        # PDBmapper accepts single or multiple protein ids
        # as input as well as prot ids stored in a file
        for ids in args.prot_id:
            # check if input is a file
            if isfile(ids) == "yes":
                input = 'file'
                with open(ids) as list_prot_ids:
                    logger.info(
                        'Input positions file contains a list of protein ids to process.')
                    # for every ensembl id
                    Parallel(n_jobs=num_cores)(delayed(wrapper)(prot_id.replace('\n', ''),
                                                                args.psdb,
                                                                args.vardb,
                                                                args.out,
                                                                args.pident,
                                                                args.evalue,
                                                                args.isoform,
                                                                args.consequence,
                                                                args.loc,
                                                                index_file,
                                                                args.dict_geneprot,
                                                               # args.uniprot,
                                                                None,
                                                                args.csv,
                                                                args.hdf)
                                               for prot_id in list_prot_ids)

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
            # print(args.prot_id)
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
                                                        args.dict_geneprot,
                                                       # args.uniprot,
                                                        None,
                                                        args.csv,
                                                        args.hdf)
                                       for ids in args.prot_id)

        # Compute execution time
        finish_message(logger, report, time_format, start, spinner)
