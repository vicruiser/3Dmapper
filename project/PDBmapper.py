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
    VEP_file : str
        Ensemble protein id 
    geneID : str
        DESCRIPTION MISSING!!
    Returns
    -------
    VCF_subset
        DESCRIPTION MISSING!!
    """
    cols = pd.read_csv(VEP_file, nrows = 0, skiprows = 42,  sep = "\t").columns
    VEP_file= open(VEP_file, 'r')
    VCF_subset = mt.parser(input_file = VEP_file,
                    ensemblID = geneID,
                    colnames = cols,
                    sep = "\t" )
    return VCF_subset

def PDBmapper(protID, interfacesDB, VCF_subset, output_dir):
    """Generate setID file.
    
    Parameters
    ----------
    protID : str
        input path of the VEP file
    interfacesDB : str
        chosen name of the output file
    VCF_subset : str
        chosen name of the output file

    Returns
    -------
    setID.File
        write data frame to a txt file with two columns. One is the gene ids and the other one the VEP file
    MappedVariants.File
        more into
    """
        # Merge them both files
    df = pd.merge(VCF_subset, interfacesDB,on=["Protein_position"],how='inner')
    setID_file = df[["region_id",
                    '#Uploaded_variation']]
    print(output_dir)
    setID_file = setID_file.drop_duplicates()
    # Save the merged dataframe
    df.to_csv(output_dir + 'MappedVariants.File', index=False, header= True,   sep = " " )
    setID_file.to_csv(output_dir + 'setID.File', index=False, header= False, sep = " " )

# define main function to execute the previous defined functions together
def main():
    # get command line options
    args = parse_argv.parse_commandline()
    
    # set interfacesDB_susbet:
    if args.intdb:
        try:
            interfacesDB = open(args.intdb, 'r')
        except IOError:
            print ("ERROR: cannot open or read input interfaces db file.")
            exit(-1)
    else: 
        print ("Default interfaces DB is used.")
        interfacesDB_file = open("./dbs/interfaces_mapped_to_v94.csv", 'r')
        interfacesDB_subset = interfaceParse(interfacesDB_file, args.protid)
    
    # set VCF_subset:
    if args.vep: 
        print("""VEP option will be available soon. Using VarMap db instead.
        Otherwise, provide your own vcf file with the -vcf option, please.""")
        args.annovar = "varmap"
    elif args.vcf:
        try: 
            print ("vcf. file provided")
            biomartdb = open("./dbs/gene_transcript_protein_ens_ids.txt", "r")
            print ("Biomart file read...")
            # get geneID
            ensemblIDs = mt.ensemblID_translator(biomartdb, args.protid)
            geneID = ensemblIDs['geneID']
            print ("We have the gene ID...", geneID)
            # get the corresponding VEP file location
            crossref_file = open("dbs/geneids_VEPfiles_crossref.txt", "r")
            print ("We have the crossref file...", crossref_file)
            # read VEP file
            VEP_filename = mt.VEP_getter(crossref_file, geneID)
            print("We have the VEP filename", VEP_filename)
            VEP_dir= args.vcf
            VEP_file = VEP_dir + '/' + VEP_filename
            print("We have the VEP file", VEP_file)
            # get subset VEP file
            VCF_subset = vcfParser(VEP_file, geneID)
        except IOError:
            print ("ERROR: cannot open or read input vcf file:")
            exit(-1)
    elif args.varmap: 
        try: 
            print ("VarMap db is used.")
            ClinVarDB =  open("./dbs/ClinVar.tsv", 'r')
            VCF_subset = vcfParser(ClinVarDB, args.protid)
        except IOError:
            print ("ERROR: cannot open or read input VarMap db file:")
            exit(-1)
    
    # set default output dir:
    if args.out is None:
        args.out = "./out/"
        
    # run PDBmapper
    PDBmapper(args.protid, interfacesDB_subset, VCF_subset, args.out)


##########################
# execute main function  #
##########################
if __name__ == '__main__':
    start = timer()
    main()
    end = timer()
    print(end - start) 



