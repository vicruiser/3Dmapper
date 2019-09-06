#!/usr/bin/python3
import argparse


def parse_commandline():

    description = \
        "%(prog)s maps rare variants to protein interfaces data in 3D."
    epilog = \
        "Copyright (c) 2019 Victoria Ruiz -- "\
        "vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel"\
        ""
    # innit parser
    parser = argparse.ArgumentParser( epilog = epilog, description = description )                 
    
    # vcf file
    annovar_group = parser.add_mutually_exclusive_group()
    annovar_group.add_argument('-vcf',   metavar = "<file>",  
                            help='input variant annotations (vcf)')
    annovar_group.add_argument("-vep",   metavar="<file>",  
                            help="input variant file to annotate them with VEP")
    annovar_group.add_argument("-varmap", action='store_true',
                            help="use VarMap db of annotated variants", )
    parser.set_defaults(annovar = "varmap")
    
    # interfaces database file
    parser.add_argument("--intdb",  dest="intdb", metavar="<file>",
                    help="Interfaces database")
    parser.set_defaults(intdb=None)

    # protein id string
    parser.add_argument("--protid",  metavar = "<String>"  ,
                    help="Ensembl protein id", required = True)
    
    # store arguments into variable
    args = parser.parse_args()

    # # clean up (recommended):
    del(parser)
    
    return args
parse_commandline()