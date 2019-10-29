#!/usr/bin/python3

# import necessary packages
import vcfpy
from scripts.decorator import tags

# add decorator to monitor function 
@tags(text_start = "Detect input annotated genomic variants file format...\n",
      text_succeed = "Detect input annotated genomic variants file format...done.\n",
      text_fail = '''Detect input annotated genomic variants file format...failed!. \n \
Please provide an file either in VCF or VEP format. \n''',
      emoji = "\U0001F50D ")

def detect_format(infile):
    '''
    Parse input and detect whether is a VCF or VEP file. Any other format
    is invalid. 
    
    Parameters
    ----------
    infile : str        
        Path to infile.
    Returns
    -------
    infile_format : str 
        Detected file format (.vep or .vcf)
    '''
    try:
        # if vcf file, check if annotated with VEP
        vcf_reader = vcfpy.Reader.from_path(infile)
        # check if "CSQ" INFO field in file
        if 'CSQ' not in vcf_reader.header.info_ids():
            # input not recognized
            vcf_reader.close()
            raise IOError()
        else: 
            # is vcf format
            return "vcf"
    except:
        # define the columns that must be present in a VEP file.  
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
        vep_file = open(infile, 'r')
        # read lines to check if header exists  
        for line in vep_file:
            res = vep_header.issubset(set(line.split()))
            # if header found, is VEP format
            if res:
                return "vep"
        # input not recognized
        raise IOError()