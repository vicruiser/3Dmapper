#!/usr/bin/env python3
# coding: utf-8

# Import necesary modules
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
import glob
from subprocess import call
import subprocess
import vcfpy
import time
from halo import Halo
import emoji

# import functions from scripts
from scripts.parse_argv import parse_commandline
from scripts.split import split 
from scripts.detect_vcf_format import detect_format
from scripts.vcf_to_vep import vcf_to_vep
from scripts.add_header import add_header
from scripts.translate_ensembl import translate_ensembl
from scripts.PDBmapper import PDBmapper
from scripts.decorator import tags

# aesthetics
description = '''

 ----------------------------------------------- Welcome to ------------------------------------------------------- 

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
        
        
------------------------  Map annotated genomic variants to protein interfaces data in 3D. ------------------------

'''


# Emojis
DNA = '\U0001F52C'  
searching_girl = '\U0001F575'

#spinner
spinner = Halo(text='Loading', spinner='dots', color= "yellow")

# define main function to execute the previous defined functions together
def main():

    # parse command line options
    args = parse_commandline()
    
    # print ascii art
    print(description)

    # create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        spinner.info(text = "Directory "+ args.out +  " created.")
    else:
        spinner.info(text = args.out + " is an existing directory. Results will be written in there.")
        
    # Manage all possible genomic variant input files
    # 1) vep = run VEP
    if args.vep:
        spinner.warn( text = 'VEP option will be available soon. Using VarMap db instead. \
        Otherwise, please provide your own vcf file with the -vcf option.')
        try:
            # ./vep -i input.vcf -o out.txt -offline
            print('VarMap db is used.')
            ClinVarDB = './default_input_data/splitted_ClinVar'
            VCF_subset = mt.parser(ClinVarDB, geneID, "\t", 'varmap')
        except:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)
    

    # 2) vcf = an annotated variant file as input either in .vcf or .vep format.
    elif args.vcf:
        # for loop in case we have multiple inputs to read
        for f in args.vcf:
            # check if file is in the right format
            try:
                print(f)
                print("LA VICKT")
                # detect the format of the vcf file(s), either .vcf or .vep
                input_format = detect_format(f)
                # set out dir and out file names
                out_dir = args.out + '/pdbmapper/input/'  # created by default
                out_file = out_dir + 'converted_vcf.vep'
                # If vcf transform into vep format and split
                if input_format == "vcf":
                    
                    # set output dir to split vep
                    vcf_db_dir = out_dir + 'vcf_db/'  # created by default

                    if os.path.isfile(out_file) is False:
                        # from vcf to vep
                        print("how many times heheheheh")
                        vcf_to_vep(f, out_dir, out_file, args.force)
                        # add header to resulting vep file
                        add_header(out_file)
                        # split vep file by protein id to speed up the
                        try:
                            split(out_file, vcf_db_dir, args.force, 'ENSG', 'vep')
                        except IOError:
                            # print error message
                            print("File " + f + " does not contain any ENSP id.")
                # If vep, only split
                elif input_format == "vep":
                    # split vep file by protein id to speed up the
                    try:
                        # set output dir to split vep
                        vcf_db_dir = './out/pdbmapper/input/vcf_db/'
                        # print starting message
                        print("Splitting VEP file. \
                        This might take a while....")
                        split(f, vcf_db_dir, args.force, 'ENSG', 'vep')
                        # print ending message
                        print("File successfully splitted.")
                    except IOError:
                        # print error message
                        print("File " + f + " does not contain any ENSP id.")

            except IOError:
                print('ERROR: cannot open or read input vcf file.')
                exit(-1)

    # 3) varmap = use VarMap as reference annotated variants file 
    elif args.varmap:
        try:
            print('VarMap db is used.')
            ClinVarDB = './dbs/splitted_ClinVar/'
            VCF_subset = mt.parser(ClinVarDB, geneID, "\t", 'varmap')

        except IOError:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)

    # set interfacesDB_susbet:
    if args.intdb:
        try:
            # set outdir
            int_db_dir = "/home/vruizser/PhD/2018-2019/git/PDBmapper/test/out/pdbmapper/input/interface_db/"
            # split interface db
            split(args.intdb, int_db_dir, args.force, 'ENSP', 'txt')
        
        except IOError:
            print('ERROR: cannot open or read input interfaces db file.')
            exit(-1)
    else:
        spinner.info(text = 'Default interfaces DB is used.')
        #set default interfaces database
        int_db_dir = "/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_interfaces_db/"

    # create chimera scripts:
    if args.chimera is not None:
        # chimera()
        pass

    # run PDBmapper
    if args.protid:
        # measure execution time
        biomartdb = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/gene_transcript_protein_ens_ids.txt'
        start = timer()
        #@tags(text_start = "Runing PDBmapper...",
        #text_succeed = "Runing PDBmapper...done.",
        #text_fail = "Runing PDBmapper...failed!",
        #emoji = "🦸")
        #def f():
        for pid in args.protid:
            # prot ids are stored in an file as list
            try:
                with open(pid) as f:
                    lines = f.read().splitlines()
                    for pids in lines:
                        # get geneID
                        try:
                            ensemblIDs = translate_ensembl(biomartdb, pids)
                            geneID = ensemblIDs['geneID']
                            try: 
                                #run PDBmapper
                                PDBmapper(pids, geneID, int_db_dir, vcf_db_dir, args.out)
                            except IOError:
                                next
                        except IOError: 
                            next
            except:
                
                # get geneID
                try:
                    print("estoy aqui")
                    print(pid)
                    print("HOLAAAAAAA!!")
                    ensemblIDs = translate_ensembl(biomartdb, pid)
                    geneID = ensemblIDs['geneID']
                    print(geneID)
                    #run PDBmapper
                    try: 
                        pass
                        #PDBmapper(pid, geneID, int_db_dir, vcf_db_dir, args.out)
                    except IOError: 
                        next
                except IOError: 
                    next
        #f()


        end = timer()
        finish = end - start
        spinner.stop_and_persist(symbol = '🧬' ,
        text = 'Congratulations!. PDBmapper has run in ' +  str(finish)  +' seconds.')
##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()

