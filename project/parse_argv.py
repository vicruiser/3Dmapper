#!/usr/bin/python3
# for commandline options:
#from optparse import OptionParser, OptionGroup
import argparse
from numpy import mean
import sys
#import urllib2
import os # This library contains functions which enables to check if a directory exists

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

    #parser = OptionParser(usage=usage, description=description,
    #                    version="%prog "+version, epilog=epilog)
    parser = argparse.ArgumentParser(prog=prog,  epilog = epilog, description = description)                 

    annovar_group = parser.add_mutually_exclusive_group(required=True)
    annovar_group.add_argument('--vcf',  dest='annovar', metavar='<file>',
                            help='input variant annotations (vcf)')
    annovar_group.add_argument("--vep",  dest="annovar", metavar="<file>",
                            help="input variant file to annotate them with VEP")
    annovar_group.add_argument("--varmap",  dest="annovar", action = "store_true",
                            help="use VarMap db of annotated variants"    )
    parser.set_defaults(annovar='varmap')

    parser.add_argument("--intdb",  dest="intdb", metavar="<file>",
                    help="Interfaces database")
    parser.set_defaults(intdb=None)

    parser.add_argument("--protid",  dest="protein_id",  action="store_true",
                    help="Ensembl protein id")
    parser.set_defaults(protein_id=False)

    args = parser.parse_args()
    # # get the options:
    if args.vcf is not None:
    # #if args.foo and args.bar is None:
    # # parser.error("--foo requires --bar. You did not specify bar.")
        print("Need an .vcf input file")
    #     # check if we have an option left (to be used as input filename):
    #     # if args:
    #     #     args.annovar = args.pop()
    #     # else:
    #     print("Need an .vcf input file")
    #     print("")
    #     parser.print_help()
    #     print("")
    #     print("ERROR: no input file given")
    #     exit(-1)

    # elif args.annovar in ('vep'):
    #     # check if we have an option left (to be used as input filename):
    #     if args:
    #         options.annovar = args.pop()
    #     else:
    #         print ("Need an input file to run VEP")
    #         print ("")
    #         parser.print_help()
    #         print ("")
    #         print ("ERROR: no input file given")
    #         exit(-1)

    # if options.intdb:
    #     # check if we have an option left (to be used as input filename):
    #     if args:
    #         options.intdb = args.pop()
    #     else:
    #         print ("Need an interfaces db input file")
    #         print ("")
    #         parser.print_help()
    #         print ("")
    #         print ("ERROR: no input file given")
    #         exit(-1)

    # # check for any leftover command line arguments:
    # if len(args):
    #     warning("ignoring additional arguments "+str(args))

    # # clean up (recommended):
    # del(parser)
    #return print(args)
    return args
    
parse_commandline()
