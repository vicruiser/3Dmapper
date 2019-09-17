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
def interfaceParse(interfacesDB, protID):
    '''Parse input interfaces database to put it in the right format.
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
    '''
    # get colnames of interfaces file
    # intDB_colnames = interfacesDB.readline().strip().split(' ')
    # parse information related to protein ID
    # prot_interface = mt.parser(input_file=interfacesDB,
    #                           ensemblID=protID,
    #    colnames=intDB_colnames,
    #    sep=' ')
    fp = './dbs/splitted_interfaces_db/' + protID + \
         '_interfaces_mapped_to_v94.csv'
    prot_interface = pd.read_csv(fp, sep=' ', header=0)

    # store subspace
    subset_prot_interface = prot_interface[['pdb.id',
                                            'ensembl.prot.id',
                                            'temp.chain',
                                            'int.chain',
                                            'interaction',
                                            'resid_sseq',
                                            'mapped.real.pos',
                                            'pdb.pos']]
    # put it into right format
    subset_prot_interface.columns = \
        subset_prot_interface.columns.str.replace('\\.', '_')

    for col in ('resid_sseq', 'mapped_real_pos', 'pdb_pos'):
        subset_prot_interface.loc[:, col] = \
            subset_prot_interface[col].str.replace('-', ',')
        subset_prot_interface.loc[:, col] = \
            subset_prot_interface[col].str.split(',')

    subset_prot_interface = mt.explode(subset_prot_interface,
                                       ['resid_sseq',
                                        'mapped_real_pos',
                                        'pdb_pos'])
    subset_prot_interface.rename(columns={'mapped_real_pos':
                                          'Protein_position'}, inplace=True)
    
    subset_prot_interface['region_id'] = subset_prot_interface['pdb_id'] + \
        '_' + subset_prot_interface['ensembl_prot_id'] + '_' + \
        subset_prot_interface['temp_chain'] + '_' + \
        subset_prot_interface['int_chain'] + '_' + \
        subset_prot_interface['interaction']
    
    return subset_prot_interface


def vcfParser(VCF_file, geneID, *args):
    '''Parse input interfaces database to put it in the right format.

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
    '''

    if args[0] == 'vcf':

        # cols = pd.read_csv(VCF_file, nrows=0, skiprows=42, sep='\t').columns
        # VCF_file = open(VCF_file, 'r')
        # VCF_subset = mt.parser(input_file=VCF_file,
        #                        ensemblID=geneID,
        #                        colnames=cols,
        #                        sep='\t')
        VCF_file = './dbs/splitted_vep_db/' + geneID + \
         '_vep_spplited.csv'
        VCF_subset = pd.read_csv(VCF_file, sep=' ', header=0)
        return VCF_subset

    elif args[0] == 'varmap':
        # cols = pd.read_csv(VCF_file, nrows=0, sep='\t').columns
        # VCF_file = open(VCF_file, 'r')
        # VCF_subset = mt.parser(input_file=VCF_file,
        #                        ensemblID=geneID,
        #                        colnames=cols,
        #                        sep='\t')
        print(geneID)
        VCF_file = './dbs/splitted_ClinVar/' + geneID + \
                   '_splitted_clinvar.csv' 
        VCF_subset = pd.read_csv(VCF_file, sep='\t', header=0)
        return VCF_subset
        # drop columns
        VCF_subset = VCF_subset[['CHROMOSOME',
                                 'COORDS',
                                 'USER_BASE',
                                 'USER_VARIANT',
                                 'ENSEMBL_BASE',
                                 'VEP_CODING_BASE',
                                 'GENE',
                                 'GENE_ACC',
                                 'TRANSCRIPT',
                                 'CODON_CHANGE',
                                 'VEP_AA',
                                 'UNIPROT_AA',
                                 'AA_CHANGE',
                                 'CHANGE_TYPE',
                                 'RES_NAME',
                                 'RES_NUM']]
        VCF_subset.drop_duplicates()
        VCF_subset = VCF_subset.rename(columns={'RES_NUM': 'Protein_position'})
        VCF_subset['#Uploaded_variation'] = VCF_subset['CHROMOSOME'] + '_' + \
            VCF_subset['COORDS'] + '_' + VCF_subset['USER_BASE'] + '_' +\
            VCF_subset['USER_VARIANT']

        return VCF_subset
    else:
        print('Error, wrong input vcf file.')


def PDBmapper(protID, interfacesDB, VCF_subset, output_dir):
    '''Generate setID file.
  
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
        write data frame to a txt file with two columns. One is the gene ids \
             and the other one the VEP file
    MappedVariants.File
        more into
    '''

    # Merge them both files
    mapped_variants = pd.concat([VCF_subset, interfacesDB],
                                axis=1, join='inner')
    # stop if there are no results
    if mapped_variants.empty:
        print('Warning:', protID, 'does not map with any annotated variant.')
        exit(-1)
    # if merging was successful, create setID file and
    # save the merged dataframe as well
    else:
        setID_file = mapped_variants[['region_id',
                                      '#Uploaded_variation']]
        setID_file = setID_file.drop_duplicates()
        # Save the merged dataframe, appending results and not
        #  reapeting headers
        with open(output_dir + 'setID.File', 'a') as f:
            setID_file.to_csv(f, sep=' ', index=False,  header=f.tell() == 0)
        with open(output_dir + 'MappedVariants.File', 'a') as f:
            mapped_variants.to_csv(f, sep=' ', index=False,
                                   header=f.tell() == 0)


# define main function to execute the previous defined functions together
def main():
    # get command line options
    args = parse_argv.parse_commandline()

    # set VCF_subset:
    if args.vep:
        print('''VEP option will be available soon. Using VarMap db instead.
        Otherwise, please provide your own vcf file with the -vcf option.''')
        try:
            print('VarMap db is used.')
            ClinVarDB = './dbs/ClinVar'
            VCF_subset = vcfParser(ClinVarDB, args.protid, 'varmap')
        except IOError:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)
    elif args.vcf:
        try: 
            print('vcf. file provided')
            biomartdb = open('./dbs/gene_transcript_protein_ens_ids.txt', 'r')
            print('Biomart file read...')
            # get geneID
            ensemblIDs = mt.ensemblID_translator(biomartdb, args.protid)
            geneID = ensemblIDs['geneID']
            print('We have the gene ID...', geneID)
            # get the corresponding VEP file location
            crossref_file = open('dbs/geneids_VEPfiles_crossref.txt', 'r')
            print('We have the crossref file...', crossref_file)
            # read VEP file
            VEP_filename = mt.VEP_getter(crossref_file, geneID)
            print('We have the VEP filename', VEP_filename)
            VEP_dir = args.vcf
            VEP_file = VEP_dir + '/' + VEP_filename
            print('We have the VEP file', VEP_file)
            # get subset VEP file
            VCF_subset = vcfParser(VEP_file, geneID, 'vcf')
        except IOError:
            print('ERROR: cannot open or read input vcf file.')
            exit(-1)
    elif args.varmap:
        try:
            print('VarMap db is used.')
            ClinVarDB = './dbs/ClinVar'
            VCF_subset = vcfParser(ClinVarDB, args.protid, 'varmap')
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
        interfacesDB_file = open('./dbs/interfaces_mapped_to_v94.csv', 'r')
        interfacesDB_subset = interfaceParse(interfacesDB_file, args.protid) 
    # set default output dir:
    if args.out is None:
        args.out = './out/'

    # create chimera scripts:
    if args.chimera is not None:
        #chimera()
        pass
    # run PDBmapper
    if args.protid:
        try:
            PDBmapper(args.protid, interfacesDB_subset, VCF_subset, args.out)
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
