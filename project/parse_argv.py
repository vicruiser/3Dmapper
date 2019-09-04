#!/usr/bin/python3
# for commandline options:
from optparse import OptionParser, OptionGroup
import argparse
from numpy import mean
import sys
import urllib2
import os # This library contains functions which enables to check if a directory exists

def parse_commandline():
    usage = "%prog <fasta> [options]"
    version = "0.0a"
    description = \
        "%prog map ."
    epilog = \
        "Copyright (c) 2019 Victoria Ruiz -- "\
        "vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel"


    parser = OptionParser(usage=usage, description=description,
                        version="%prog "+version, epilog=epilog)

    mutex_group = parser.add_mutually_exclusive_group()
    parser.add_option("", "--vcf",  dest="vcf", metavar="<file>",
                    help="input variant annotations (vcf)")
    parser.set_defaults(vcf=None)
    
    parser.add_option("", "--vep",  dest="vep", metavar="<file>",
                    help="input variant file to annotate them with VEP")
    parser.set_defaults(vep=None)

    parser.add_option("-vm", "--varmap",  dest="varmap", action = "store_true",
                    help="use VarMap db of annotated variants")
    parser.set_defaults(varmap=True)



mutex_group.add_argument("--show", action="store_const", 
    dest="mutex", const="show")
mutex_group.add_argument("--insert", action="store_const", 
    dest="mutex", const="insert")
mutex_group.add_argument('--delete', action="store_const",
    dest="mutex", const="delete")


parser.set_defaults(mutex='show')
args = parser.parse_args()
print(args)


    parser.add_option("", "--intdb",  dest="intdb", metavar="<file>",
                    help="Interfaces database")
    parser.set_defaults(int_db=None)

    parser.add_option("-p", "--protid",  dest="protein_id",  action="store_true",
                    help="Ensembl protein id")
    parser.set_defaults(align_local=False)


    # get the options:
    (options, args) = parser.parse_args()

    if not options.vcf:
        # check if we have an option left (to be used as input filename):
        if args:
            options.varmap = args.pop()
        elif args: 
            options.vep = args.pop()
        else:
            print ("Need at least an input file (--vcf , --varmap or vcf)")
            print ("")
            parser.print_help()
            print ("")
            print ("ERROR: no input file given")
            exit(-1)

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))

    # clean up (recommended):
    del(parser)
    return options
    