#!/usr/bin/python3
import argparse, os

def dir_path(dirName):
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
    else:    
        print("Directory " , dirName ,  " already exists")

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
    annovar_group.add_argument('-vcf',   metavar = "<String>",  dest = "vcf",
                            help='input directory containing all annotated variants files (vcf)')
    annovar_group.add_argument("-vep",   metavar="<file>",  dest = "vep",
                            help="default VEP input")
    annovar_group.add_argument("-varmap", action='store_true', dest = "varmap",
                            help="use ClinVar db of annotated variants", default = "varmap" )
    #parser.set_defaults(annovar = "varmap")
    
    # interfaces database file
    parser.add_argument("-intdb",  dest="intdb", metavar="<file>",
                    help="Interfaces database")
    parser.set_defaults(intdb=None)

    # protein id string
    parser.add_argument("-protid",  metavar = "<String>"  , dest= "protid",
                    help="Ensembl protein id", required = True)
    
    # create default output directory
    parser.add_argument("-out",  metavar = "<String>"  , dest = "out",
                    help="output directory", type=dir_path )
    parser.set_defaults(out="./out/")
    # store arguments into variable
    args = parser.parse_args()

    # # clean up (recommended):
    del(parser)
    
    return args





