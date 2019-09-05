#!/usr/bin/python3
# for commandline options:
#from optparse import OptionParser, OptionGroup
import argparse
from numpy import mean
import sys
#import urllib2
import os # This library contains functions which enables to check if a directory exists
import re
def parse_commandline():
    
    prog='pdb_mapper'
    usage = "%prog <fasta> [options]"
    version = "0.0a"
    description = \
        "%prog map ."
    epilog = \
        "Copyright (c) 2019 Victoria Ruiz -- "\
        "vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel"\
        ""

    parser = argparse.ArgumentParser(prog=prog,  epilog = epilog, description = description)                 

    annovar_group = parser.add_mutually_exclusive_group()
    annovar_group.add_argument('--vcf',  dest = "annovar",  
                            help='input variant annotations (vcf)')
    annovar_group.add_argument("--vep",   metavar="<file>", dest = "annovar",  
                            help="input variant file to annotate them with VEP")
    annovar_group.add_argument("--varmap", dest = "annovar", default = 'annovar',
                            help="use VarMap db of annotated variants", nargs='?',
                            const='value_if_noarg')
    parser.set_defaults(annovar = "varmap")

    parser.add_argument("--intdb",  dest="intdb", metavar="<file>",
                    help="Interfaces database")
    parser.set_defaults(intdb=None)

    parser.add_argument("--protid",  dest="protein_id",  action="store_true",
                    help="Ensembl protein id")
    parser.set_defaults(protein_id=False)
    
    args = parser.parse_args()

    # raise an exception if we weren't able to find a match
    if args.intdb is None:
        print ("Default interfaces DB is used.")

    if args.annovar == "varmap":
        print ("VarMap db is used.")

    print(args.annovar)

    # # clean up (recommended):
    del(parser)
    #return print(args)
    return args
