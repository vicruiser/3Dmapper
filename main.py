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
from scripts.translate_ensembl import translate_ensembl
from scripts.PDBmapper import PDBmapper

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
                # detect the format of the input file
                print('Detecting inputs variants file format...')
                # detect the format of the vcf file(s), either .vcf or .vep
                input_format = detect_format(f)
                # set out dir and out file names
                out_dir = args.out + '/pdbmapper/input/'  # created by default
                out_file = out_dir + 'converted_vcf.vep'

                # If vcf transform into vep format and split
                if input_format == "vcf":

                    # print starting message
                    print("The file is in .vcf format. Transforming into VEP format. \
                        This might take up some time....")
                    # from vcf to vep
                    vcf_to_vep(f, out_dir, out_file, args.force)
                    # print ending message
                    print("File succesfully converted.")

                    # add header to resulting vep file
                    add_header = "gawk -i inplace '{{if(NR==1){{$0=\"#Uploaded_variation Location Allele Gene \
    Feature Feature_type Consequence cDNA_position CDS_position Protein_position \
    Amino_acids Codons Existing_variation\\n\"$0; print $0}}; if(NR!=1){{print $0}}}}' {}"
                    cmd2 = add_header.format(out_file)
                    p2 = subprocess.Popen(
                        cmd2, stdout=subprocess.PIPE, shell=True)

                    # set output dir to split vep
                    vcf_db_dir = out_dir + 'vcf_db/'  # created by default
                    # print starting message
                    print("Splitting VEP file. \
                        This might take a while....")
                    # split vep file by protein id to speed up the
                    try:
                        split(out_file, vcf_db_dir, args.force, 'ENSG', 'vep')
                        # print ending message
                        print("File successfully splitted.")
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
            # print starting message
            print("Splitting interfaces file...")
            # split interface db
            split(args.intdb, int_db_dir, args.force, 'ENSP', 'txt')
            
            print("finished!")
        except IOError:
            print('ERROR: cannot open or read input interfaces db file.')
            exit(-1)
    else:

        spinner.start('Default interfaces DB is used.')
        
        #set default interfaces database
        int_db_dir = "/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/splitted_interfaces_db/"
                      
        spinner.stop_and_persist(symbol= DNA + searching_girl ,text= ' Default interfaces DB is used.')

    # create chimera scripts:
    if args.chimera is not None:
        # chimera()
        pass

    # run PDBmapper
    if args.protid:
        # measure execution time
        spinner.start("Running PDBmapper")
        start = timer()
        print(args.protid)
        
        # prot ids are stored in an file as list
        try:
            with open(args.protid) as f:
                lines = f.read()
                for line in lines:
                    # get geneID
                    biomartdb = open('/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/gene_transcript_protein_ens_ids.txt', 'r')
                    spinner.start('Reading Biomart...')
                    ensemblIDs = translate_ensembl(biomartdb, pid)
                    geneID = ensemblIDs['geneID']
                    spinner.succeed('ENSP translated to ENSG')
                    print(geneID)
                    print(pid, geneID, int_db_dir, vcf_db_dir, args.out)
                    PDBmapper(pid, geneID, int_db_dir, vcf_db_dir, args.out)
                # or

        except:
            for pid in args.protid:
                # get geneID
                biomartdb = open('/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/gene_transcript_protein_ens_ids.txt', 'r')
                spinner.start('Reading Biomart...')
                ensemblIDs = translate_ensembl(biomartdb, pid)
                geneID = ensemblIDs['geneID']
                spinner.succeed('ENSP translated to ENSG')
                print(geneID)


                print(pid, geneID, int_db_dir, vcf_db_dir, args.out)
                PDBmapper(pid, geneID, int_db_dir, vcf_db_dir, args.out)

        #                 except IOError:
        #                     print('ERROR: Ensembl protein id provided is no supported.')
        #                     exit(-1)

            finish = end - start
            #print('Congratulations!. PDBmapper has run in', finish, 'seconds.')
            spinner.succeed(text = 'Congratulations!. PDBmapper has run in  seconds.')
##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()

