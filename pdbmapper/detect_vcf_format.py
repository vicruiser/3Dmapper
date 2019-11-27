# coding: utf-8

# import necessary packages
import vcfpy
from .decorator import tags

# add decorator to monitor function


@tags(text_start="Detect input annotated genomic variants file format...\n",
      text_succeed="Detect input annotated genomic variants file format...done.\n",
      text_fail='''Detect input annotated genomic variants file format...failed!. \n \
Please provide a file either in VCF or VEP format. \n''',
      emoji="\U0001F50D ")
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
    str 
        Detected file format ("vep" or "vcf")
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
        # define the columns that must be present in the VEP file.
        vep_header = set(['#Uploaded_variation',
                          'Gene',
                          'Feature',
                          'Consequence',
                          'Protein_position'
                          ])
        # define the columns that must be present in an alternative variant file.
        alt_header = set(['Uploaded_variation',
                          'Gene',
                          'Feature',
                          'Consequence',
                          'Protein_position'
                          ])
        # open file to parse it
        variants_file = open(infile, 'r')
        # read lines to check if header exists
        for line in variants_file:
            isvep = vep_header.issubset(set(line.split()))
            # if header found, is VEP format
            if isvep:
                return "vep"
            # if vep format-like
            isalt = alt_header.issubset(set(line.split()))
            if isalt:
                return "alt"
        # input not recognized
        raise IOError()
