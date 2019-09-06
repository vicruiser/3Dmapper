#!/usr/bin/env python3
# coding: utf-8

# Import necesary modules 
import mapping_tools as mt
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
from VEPcrossref import VEPfileCrossrefGenerator as cr
import parse_argv

# Extract the info corresponding to the prot ID (Interface parse)

#protID = 'ENSP00000482258'
#int_db_file = './dbs/interfaces_mapped_to_v94.csv'

def interfaceParse (interfacesDB, protID):
    """Parse input interfaces database to put it in the right format.
    
    Parameters
    ----------
    protID : str
        Ensemble protein id 
    interfacesDB_filepath : str
        DESCRIPTION MISSING!!
    Returns
    -------
    subset_interfaces_db
        DESCRIPTION MISSING!!
    """
    # get colnames of interfaces file
    intDB_colnames = interfacesDB.readline().strip().split(" ")
    # parse information related to protein ID
    prot_interface = mt.parser(input_file = interfacesDB,
                                ensemblID = protID,
                                colnames =  intDB_colnames,
                                sep = " " )
    # store 
    subset_prot_interface = prot_interface[['pdb.id',
                                            'ensembl.prot.id',
                                            'temp.chain',
                                            'int.chain',
                                            'interaction',
                                            'resid_sseq',
                                            'mapped.real.pos',
                                            'pdb.pos' ]]
    # put it into right format
    subset_prot_interface.columns = subset_prot_interface.columns.str.replace("\\.", "_")

    for col in ("resid_sseq", "mapped_real_pos", "pdb_pos"):
        subset_prot_interface.loc[ : , col] = subset_prot_interface[col].str.replace("-", ',')
        subset_prot_interface.loc[ : , col] = subset_prot_interface[col].str.split(',')

    subset_prot_interface = mt.explode(subset_prot_interface, ["resid_sseq", "mapped_real_pos", "pdb_pos"])
    subset_prot_interface.rename(columns={'mapped_real_pos':'Protein_position'}, inplace=True)
    subset_prot_interface['region_id'] = subset_prot_interface['pdb_id'] + '_' + subset_prot_interface['ensembl_prot_id'] + '_' + subset_prot_interface['temp_chain'] + '_' + subset_prot_interface['int_chain'] + '_' + subset_prot_interface['interaction']
    
    return subset_prot_interface

def vcfParser(VEP_file, geneID):
    """Parse input interfaces database to put it in the right format.
    
    Parameters
    ----------
    protID : str
        Ensemble protein id 
    interfacesDB_filepath : str
        DESCRIPTION MISSING!!
    Returns
    -------
    subset_interfaces_db
        DESCRIPTION MISSING!!
    """

    
    VEP_dir= "/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/"
    geneID = ensemblIDs['geneID']
    VEPfile = open(VEP_dir + VEP_filename, "r")


    cols = pd.read_csv(VEPfile, nrows = 0, skiprows = 42, sep = "\t").columns
    VCF = mt.parser(input_file = VEPfile,
                    ensemblID = geneID,
                    colnames = cols,
                    sep = "\t" )

    return VCF

def PDBmapper(protID, interfacesDB, annovar):
    """Generate setID file.
    
    Parameters
    ----------
    protID : str
        input path of the VEP file
    interfacesDB_filepath : str
        chosen name of the output file

    Returns
    -------
    setID.File
        write data frame to a txt file with two columns: one is the gene ids and the other one the VEP file
    data.frame 
        kasklkssklasdkmaskldmakmdalkmd
    """

    setID = mt.MapVariantToPDB(annovar, interfacesDB, 'region_id')
    setID.to_csv("setID.File")

# define main function to execute the previous defined functions together
def main():

    # get command line options
    args = parse_argv.parse_commandline()
    
    # set interfaces db:
    if args.intdb:
        try:
            interfacesDB = open(args.intdb, 'r')
        except IOError:
            print ("ERROR: cannot open or read input interfaces db file.")
            exit(-1)
    else: 
        print ("Default interfaces DB is used.")
        interfacesDB_file = open("./dbs/interfaces_mapped_to_v94.csv", 'r')
        interfacesDB = interfaceParse(interfacesDB_file, args.protid)
    
    # set annovar:
    if args.annovar == "vep": 
        print("""VEP option will be available soon. Using VarMap db instead.
        Otherwise, provide your own vcf file with the -vcf option, please.""")
        annovar_file =  open("output_vep.txt", 'r')
    elif args.annovar == "vcf":
        try: 
            print ("VarMap db is used.")
            biomartdb = open("./dbs/gene_transcript_protein_ens_ids.txt", "r")
            ensemblIDs = mt.ensemblID_translator(biomartdb, args.protID)
            # get the corresponding VEP file
            crossref_file = open("dbs/gene_transcript_protein_ens_ids.txt", "r")
            VEP_filename = mt.VEP_getter(crossref_file, ensemblIDs["geneID"])
            annovar_file = open(args.annovar, 'r')
            
            annovar = vcfParser(annovar_file, args.protID)
        except IOError:
            print ("ERROR: cannot open or read input vcf file:")
            exit(-1)
    elif args.annovar == "varmap": 
        try: 
            print ("VarMap db is used.")
            annovar_file =  open("./PDBmapper/dbs/VarMap.csv", 'r')
            annovar = vcfParser(annovar_file, args.protID)
        except IOError:
            print ("ERROR: cannot open or read input VarMap db file:")
            exit(-1)
    
    # run PDBmapper
    #PDBmapper(args.protid, interfacesDB, annovar)



##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()
