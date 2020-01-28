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
import multiprocessing as mp

from halo import Halo
from timeit import default_timer as timer
from subprocess import call
from tabulate import tabulate

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


class MapTools:

    time = '[' + time.ctime(time.time()) + '] '

    def log(self, message, report, logger):

        logger.info(message)
        report.write(self.time + message + '\n')

    def run():
        pass

    def run_in_parallel():
        pass

    def stats(self, var_statsfile, int_statsfile, mapped_infofile, out_dir):
        # PDBmapper stats
        var_stats = pd.read_csv(var_statsfile, sep=" |\t", engine='python')
        int_stats = pd.read_csv(int_statsfile, sep=" |\t", engine='python')

        # number of mapped variants
        n_variants_unique_cmd = "awk 'NR>1{{print $1}}' {} | uniq | wc -l"
        n_variants_unique, err1 = call_subprocess(
            n_variants_unique_cmd.format(mapped_infofile))

        n_variants_cmd = "awk 'NR>1{{print $1}}' {}  | wc -l"
        n_variants, err1 = call_subprocess(
            n_variants_cmd.format(mapped_infofile))

        # number of mapped genes
        n_genes_cmd = "awk 'NR>1{{print $4}}' {} | uniq | wc -l"
        n_genes, err2 = call_subprocess(n_genes_cmd.format(mapped_infofile))

        # process results
        n_variants, n_variants_unique, n_genes = n_variants.decode(
            'utf-8').rstrip(), n_variants_unique.decode(
            'utf-8').rstrip(), n_genes.decode('utf-8').rstrip()

        # number of mapped interfaces, type of interfaces and consequences
        interfaces_info = pd.read_csv(mapped_infofile, sep=' ', usecols=[
                                      'interface_id', 'Consequence'])
        interfaces_info[['pdb_id', 'ENSP', 'temp_chain', 'int_chain', 'type']
                        ] = interfaces_info.interface_id.str.split("_", expand=True)

        # Consequence table counts
        summary = interfaces_info.groupby(
            ['type', 'Consequence']).size().reset_index(name='Count')
        summary['type'] = summary['type'].replace({'ligand': 'interface_with_ligand',
                                                   'protein': 'interface_with_protein',
                                                   'nucleic': 'interface_with_dna'})

        summary = summary.pivot(index='Consequence',
                                columns='type', values='Count')
        summary = summary.fillna(0)

        summary.to_csv(os.path.join(out_dir, 'pdbmapper_consequences_stats.info'),
                       sep=' ', encoding='utf-8', index=True)
        stats_table = tabulate(
            summary, headers='keys', tablefmt='psql')

        # number of mapped proteins
        n_prot = interfaces_info['ENSP'].drop_duplicates().count()

        # number of mapped interfaces
        n_int = interfaces_info['interface_id'].drop_duplicates().count()

        # number of interfaces ligand
        n_int_ligand = interfaces_info.loc[
            interfaces_info['type'] == "ligand", 'interface_id'].drop_duplicates().count()

        # number of interfaces protein
        n_int_prot = interfaces_info.loc[
            interfaces_info['type'] == "protein", 'interface_id'].drop_duplicates().count()

        # number of interfaces nucleic
        n_int_dna = interfaces_info.loc[
            interfaces_info['type'] == "nucleic", 'interface_id'].drop_duplicates().count()

        # stats table
        stats_message = {'Variable': ['Variants ids', 'Unique variants ids',
                                      'Interface ids', 'Interface ids (dna)',
                                      'Interface ids (ligand)', 'Interface ids (protein)',
                                      'Protein ids', 'Gene ids'],
                         'Total':
                         [int(n_variants), int(n_variants_unique),
                          n_int, n_int_dna,
                          n_int_ligand, n_int_prot,
                          n_prot, int(n_genes)],
                         '% (mapped / input)':
                         ['-',
                          round((int(n_variants_unique) /
                                 var_stats.iloc[0]['n_variants']) * 100, 2),
                          round((int(n_int) /
                                 int_stats.iloc[0]['n_interfaces']) * 100, 2),
                          round((int(n_int_dna) /
                                 int_stats.iloc[0]['n_interfaces_nucleic']) * 100, 2),
                          round((int(n_int_ligand) /
                                 int_stats.iloc[0]['n_interfaces_ligand']) * 100, 2),
                          round((int(n_int_prot) /
                                 int_stats.iloc[0]['n_interfaces_protein']) * 100, 2),
                          round((int(n_prot) /
                                 int_stats.iloc[0]['n_ENSP']) * 100, 2),
                          round((int(n_genes) /
                                 var_stats.iloc[0]['n_genes']) * 100, 2)]}
        # conver it to data frame and save it
        stats_message = pd.DataFrame(stats_message)
        stats_message.set_index('Variable')
        stats_message2 = tabulate(
            stats_message, headers='keys', tablefmt='psql')
        stats_message.to_csv(os.path.join(out_dir, 'pdbmapper_stats.info'),
                             sep='\t', encoding='utf-8', index=False)

        return stats_message2, stats_table


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
          -------------------------------------------------------------------------
         |  Copyright (c) 2019 Victoria Ruiz --                                    |
         |  vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          -------------------------------------------------------------------------

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
        p = mp.Pool(mp.cpu_count()-1)

    else:
        p = mp.Pool(1)
        # results = pool.map(howmany_within_range_rowonly, [row for row in data])
        # pool.close()

    if args.force is True and os.listdir(args.out):
        pass
        # for item in os.listdir(args.out):
        #    if item.endswith(".File"):
        #        os.remove(os.path.join(args.out, item))
    else:
        logger.warning(
            'Directory ' + args.out + 'is not empty. Not overwritting files. ' +
            'Please select option --force or specify a different output dir.')
        spinner.warn(
            text=' Directory ' + args.out + 'is not empty. Not overwritting files. ' +
            'Please select option --force or specify a different output dir.')
        exit(-1)

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
                try:
                    wrapper(ensemblid,
                            args.intdb,
                            args.vardb,
                            args.out,
                            args.pident,
                            args.consequence)
                except:
                    continue
            elif isfile(ids) == "not_recognized":
                pdbmapper.log('The input is neither an id(s) or a file containing a list of ids.',
                              report, logger)
                spinner.fail(
                    'The input is neither an id(s) or a file containing a list of ids.')
                exit(-1)

        # compute statistics

        var_statsfile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
            args.vardb, 'makevariantsdb_stats.info'))[0]))

        int_statsfile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
            args.intdb, 'makeinterfacesdb_stats.info'))[0]))

        mapped_infofile = os.path.abspath(os.path.normpath(glob.glob(os.path.join(
            args.out, 'MappedVariants*' + str(args.pident) + '*File'))[0]))

        stats_message, stats_table = maptools.stats(
            var_statsfile, int_statsfile, mapped_infofile, args.out)

        # Compute execution time
        end = time.time()
        logger.info('Done.')
        report.write(time_format + 'PDBmapper has finished successfully. Total time: ' +
                     str(datetime.timedelta(
                         seconds=round(end-start))) + '\n \n')
        report.write('Stats\n')
        report.write('-----\n')
        report.write(stats_message)
        report.write('\n')
        report.write(stats_table)
        report.write('\n')
        # print in console result
        spinner.stop_and_persist(symbol='\U0001F4CD',
                                 text=' PDBmapper has finished successfully. Total time: ' +
                                 str(datetime.timedelta(
                                     seconds=round(end-start))))
