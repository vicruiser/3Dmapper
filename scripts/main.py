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

# import functions from scripts
import parse_argv
from VEPcrossref import VEPfileCrossrefGenerator as cr
import mapping_tools as mt
import PDBmapper
import scripts.vep_parser
from scripts.detect_vcf_format import detect_format
from scripts.vcf_to_vep import vcf_to_vep
from scripts.split_vep import split_vep

# define main function to execute the previous defined functions together
def main():
    # get command line options
    args = parse_argv.parse_commandline()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        print("Directory", args.out,  "created.")
    else:    
        print(args.out,
              "is an existing directory. Results will be written in there.")
    # get geneID
        # biomartdb = open('./dbs/gene_transcript_protein_ens_ids.txt', 'r')
        # print('Biomart file read...')
        # ensemblIDs = mt.ensemblID_translator(biomartdb, args.protid)
        # geneID = ensemblIDs['geneID']

    # set VCF_subset:
    if args.vep:
        print('''VEP option will be available soon. Using VarMap db instead.
        Otherwise, please provide your own vcf file with the -vcf option.''')
        try:
            # ./vep -i input.vcf -o out.txt -offline
            print('VarMap db is used.')
            ClinVarDB = './dbs/splitted_ClinVar'
            VCF_subset = vcfParser(ClinVarDB, geneID, "\t", 'varmap')
        
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    elif args.vcf:
        try:
            print('Detecting inputs variants file format...')
            # for loop in case we have multiple inputs to read
            for f in args.vcf: 
                # detect the format of the vcf file(s)
                # possible formats: .vcf or .vep
                input_format = detect_format(f)
                if input_format == "vcf":
                    # 1) put vcf into vep format
                    out_dir = './pdbmapper/input/' #created by default
                    out_file = out_dir + 'converted_vcf.vep'
                    vcf_to_vep(f, out_dir, out_file)
                    # 2) split vep file by protein id to speed up the
                    # mapping proccess
                    vcf_db_dir = out_dir + 'vcf_db/' #created by default
                    split_vep(f, vcf_db_dir)
                elif input_format == "vep":
                    # split vep file by protein id to speed up the
                    # mapping proccess
                    vcf_db_dir = './pdbmapper/input/vcf_db/'
                    split_vep(f, vcf_db_dir) 
                else: 
                    print(input_format)
        except IOError:
            print('ERROR: cannot open or read input vcf file.')
            exit(-1)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    elif args.varmap:
        try:
            print('VarMap db is used.')
            ClinVarDB = './dbs/splitted_ClinVar/'
            VCF_subset = vcfParser(ClinVarDB, geneID, "\t", 'varmap')
        except IOError:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)

    # set interfacesDB_susbet:
    if args.intdb:
        try:
            interfacesDB = open(args.intdb, 'r')
        except IOError:
            print('ERROR: cannot open or read input interfaces db file.')
            exit(-1)
    else:
        print('Default interfaces DB is used.')
        interfaces_dir = './dbs/splitted_interfaces_db'
        interfacesDB_subset = interfaceParse(interfaces_dir, args.protid)


    # create chimera scripts:
    if args.chimera is not None:
        #chimera()
        pass

    # run PDBmapper
    if args.protid:
        # input list of proteins
        if isinstance(args.protid, list):
            # iterate over list
            for one_protid in args.protid:
                try:
                    PDBmapper(one_protid, interfacesDB_subset, VCF_subset,
                            args.out)
                except IOError:
                    print('ERROR: Ensembl protein id provided is no supported.')
                    exit(-1)
        # single protein id as input
        else:
            try:
                PDBmapper(args.protid, interfacesDB_subset, VCF_subset,
                          args.out)
            except IOError:
                print('ERROR: Ensembl protein id provided is no supported.')
                exit(-1)
    if args.protid_file: 
         
        for one_protid in args.protid:
            try:
                PDBmapper(one_protid, interfacesDB_subset, VCF_subset,
                        args.out)
            except IOError:
                print('ERROR: Ensembl protein id provided is no supported.')
                exit(-1)

##########################
# execute main function  #
##########################
if __name__ == '__main__':
    # measure execution time
    start = timer()
    main()
    end = timer()
    finish = end - start
    print('Congratulations!. PDBmapper has run in', finish, 'seconds.')
