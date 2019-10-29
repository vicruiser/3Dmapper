#!/usr/bin/python3

# import necessary packages
import vcfpy
from scripts.decorator import tags

@tags(text_start = "Detect inputs variants file format...\n",
      text_succeed = "Detect inputs variants file format...done.\n",
      text_fail = "Detect inputs variants file format...failed!\n",
      emoji = "\U0001F50D ")
def detect_format(infile):
    '''Parse input and detect whether is a VCF or VEP file. Any other format
    is invalid. 
    Parameters
    ----------
    infile : str        
        Path to infile.
    Returns
    -------
    infile_format : str 
        Detected file format.
    '''
    try:
        # if vcf file, check if annotated with VEP
        vcf_reader = vcfpy.Reader.from_path(infile)
        if 'CSQ' not in vcf_reader.header.info_ids():
            vcf_reader.close()
            return print("ERROR: VCF is not VEP-annotated. Please annotate\
                          the VCF with VEP before running this tool.")
        else : 
            print('Input file format is VCF')
            infile_format = "vcf"
            return infile_format
    except:
        print('The input file is not in vcf format.\
               Checking if it is a VEP file...')
        # define the columns that must be present in the vep file.  
        vep_header = set(['#Uploaded_variation',
                                'Location',
                                'Allele',
                                'Gene',
                                'Feature',
                                'Feature_type',
                                'Consequence',
                                'cDNA_position',
                                'CDS_position',
                                'Protein_position',
                                'Amino_acids',
                                'Codons',
                                'Existing_variation'])
        # open file to parse it     
        vep_file = open(infile,'r')
        # read lines to check if header exists in it 
        for line in vep_file:
            res = vep_header.issubset(set(line.split()))
            if res:
                print("is vep format") 
                infile_format = "vep"
                return infile_format

        raise IOError("Sorry, input format not recognized. Possible input\
                     formats are .vep or .vcf")