#!/usr/bin/env python3
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
from scripts.parse_argv import parse_commandline
from scripts.run_vep import run_vep
from scripts.split import split
from scripts.detect_vcf_format import detect_format
from scripts.vcf_to_vep import vcf_to_vep
from scripts.add_header import add_header
from scripts.translate_ensembl import translate_ensembl
from scripts.PDBmapper import PDBmapper
from scripts.decorator import tags

pd.options.mode.chained_assignment = None

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


# Emojis
DNA = '\U0001F9EC'
searching_girl = '\U0001F575'

# spinner
spinner = Halo(text='Loading', spinner='dots12', color="red")

# define main function


def main():

    # parse command line options
    args = parse_commandline()

    # print ascii art
    print(description)

    # create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        spinner.info(text="Directory " + args.out + " created.\n")
    else:
        spinner.info(
            text=args.out + " is an existing directory. Results will be written in there.\n")

    # Manage all possible genomic variant input files
    ####################
    # 1) vep = run VEP #
    ####################
    if args.vep:
        spinner.warn(text='VEP option will be available soon. Using VarMap db instead. \
Otherwise, please provide your own vcf file with the -vcf option.\n')
        try:
            run_vep()
        except IOError:
            vcf_db_dir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'
            spinner.info('Using VarMap db\n')
            exit(-1)
    ##############################################################################
    # 2) vcf = an annotated variant file as input either in .vcf or .vep format. #
    ##############################################################################
    elif args.vcf:
         # compute total time of splitting vcf file
        start = timer()

        # set out dir and out file names
        # created by default
        out_dir = os.path.join(args.out, 'input')
        out_file = os.path.join(
            out_dir, 'converted_vcf.vep')  # created by default
        # set output dir to split vep
        vcf_db_dir = os.path.join(out_dir, 'vcf_db')  # created by default
        # create output dir if it doesn't exist
        os.makedirs(vcf_db_dir, exist_ok=True)
        # change input format if file doesn't exists or overwrite is True
        if not os.listdir(vcf_db_dir) or args.force.lower() == 'y':
            # for loop in case we have multiple inputs to read
            for f in args.vcf:
                # check if input is a file
                try:
                    with open(f) as list_var_files:
                        var_f = list_var_files.read().splitlines()
                        # for every prot id
                        for var_path in var_f:
                            # detect the format of the vcf file(s), either .vcf or .vep
                            input_format = detect_format(var_path)
                            # If vcf transform into vep format and split
                            if input_format == "vcf":
                                # change input format if file doesn't exists or overwrite is True
                                if os.path.isfile(out_file) is False or args.force.lower() == 'y':
                                    # from vcf to vep
                                    vcf_to_vep(var_path, out_dir,
                                               out_file, args.force)
                                    # add header to resulting vep file
                                    add_header(out_file)
                                    # split vep file by protein id to speed up the
                                    # mapping process
                                    split('ENSG', out_file, vcf_db_dir,
                                          'vep', args.force)

                            # If vep, only split
                            elif input_format == "vep" or input_format == "alt":
                                # split if empty dir or overwrite is True
                                if not os.listdir(vcf_db_dir) or args.force.lower() == 'y':
                                    # split vep file by protein id to speed up the
                                    # mapping process
                                    split('ENSG', var_path, vcf_db_dir,
                                          'vep', args.force)
                            else:
                                print('Warning: input file ' + var_path +
                                      ' is not in vep nor vcf format.')
                                continue
                except:
                    # change input format if file doesn't exists or overwrite is True
                    if not os.listdir(vcf_db_dir) or args.force.lower() == 'y':
                        # detect the format of the vcf file(s), either .vcf or .vep
                        input_format = detect_format(f)
                        # If vcf transform into vep format and split
                        if input_format == "vcf":
                            # change input format if file doesn't exists or overwrite is True
                            if os.path.isfile(out_file) is False or args.force.lower() == 'y':
                                # from vcf to vep
                                vcf_to_vep(f, out_dir, out_file, args.force)
                                # add header to resulting vep file
                                add_header(out_file)
                                # split vep file by protein id to speed up the
                                # mapping process
                                split('ENSG', out_file, vcf_db_dir,
                                      'vep', args.force)

                        # If vep, only split
                        elif input_format == "vep" or input_format == "alt":
                            # split if empty dir or overwrite is True
                            if not os.listdir(vcf_db_dir) or args.force.lower() == 'y':
                                # split vep file by protein id to speed up the
                                # mapping process
                                split('ENSG', f, vcf_db_dir, 'vep', args.force)
                        else:
                            print('Warning: input file ' + var_path +
                                  ' is not in vep nor vcf format.')
                            continue
        # time execution
        end = timer()
        finish = end - start
        log_finish = open(os.path.join(out_dir, 'results_report.txt'), 'a')
        log_finish.write('The conversion and splitting of the vcf file has taken ' +
                         str(finish/60) + ' minutes.')
    ###################################################################
    # 3) maf = input file is in Mutation Annotated File (.maf) format #
    ###################################################################
    elif args.maf:

        pass
    ###############################################################
    # 4) varmap = use VarMap as reference annotated variants file #
    ###############################################################
    elif args.varmap:
        spinner.info('Using VarMap db')
        vcf_db_dir = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_ClinVar'

    # set interfacesDB_susbet:
    if args.intdb:
        # set outdir
        int_db_dir = "/home/vruizser/PhD/2018-2019/git/PDBmapper/test/out/pdbmapper/input/interface_db"
        # split interface db
        split('ENSP', args.intdb, int_db_dir, 'txt', args.force)
        # set origin of input interface
        input_intdb = 'external'
    else:
        spinner.info(text='Default interfaces DB is used.')
        # set default interfaces database
        int_db_dir = "/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_interfaces_db"
        # set origin of input interface
        input_intdb = 'default'
    # create chimera scripts:
    if args.chimera is not None:
        # chimera()
        pass

    # run PDBmapper
    if args.protid:
        # output dir
        out_dir = os.path.join(args.out, 'results')
        # create output dir if it doesn't exist
        os.makedirs(out_dir, exist_ok=True)
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
            for prot_ids in args.protid:
                # check if input is a file
                try:
                    with open(prot_ids) as f:
                        lines = f.read().splitlines()
                        # set variable input as file
                        input = "file"

                except:
                    # set variable input as not file
                    input = "not_file"

                if input == "file":
                    for prot_id in lines:
                        try:
                            # for pids in lines:
                            ensemblIDs = translate_ensembl(
                                prot_id, args.filter_iso)
                            geneID = ensemblIDs['geneID']
                        except IOError:
                            log = open(os.path.join(
                                out_dir, 'log_ensembl.File'), 'a')
                            log.write('Warning: ' + prot_id +
                                      ' has no ENGS.\n')
                            continue
                        # run PDBmapper
                        try:
                            PDBmapper(prot_id,
                                      geneID,
                                      int_db_dir,
                                      input_intdb,
                                      vcf_db_dir,
                                      out_dir,
                                      args.pident,
                                      args.filter_var)
                        # error handling
                        except IOError:
                            continue

                # input is not a file but one or more protein ids
                # given in command line
                elif input == "not_file":
                    # for prot id get the gene id
                    prot_id = prot_ids
                    try:
                        ensemblIDs = translate_ensembl(
                            prot_id, args.filter_iso)
                        geneID = ensemblIDs['geneID']
                        # run PDBmapper
                        try:
                            PDBmapper(prot_id,
                                      geneID,
                                      int_db_dir,
                                      input_intdb,
                                      vcf_db_dir,
                                      out_dir,
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

        log_finish = open(os.path.join(out_dir, 'results_report.txt'), 'a')
        log_finish.write('Congratulations!. PDBmapper has run in ' +
                         str(finish) + ' minutes.')


##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()
