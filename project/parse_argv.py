#!/usr/bin/python3
# for commandline options:
from optparse import OptionParser, OptionGroup
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

    # sequence/alignment options:
    parser.add_option("-f", "--fasta",  dest="fasta", metavar="<file>",
                    help="input alignment file (fasta)")
    parser.set_defaults(fasta=None)

    parser.add_option("-e", "",  dest="exchange_matrix",
                    help="Exchange matrix: pam250, blosum62 or identity (%default)")
    parser.set_defaults(exchange_matrix="pam250")

    parser.add_option("-l", "",  dest="align_local",  action="store_true",
                    help="align local")
    parser.set_defaults(align_local=False)

    parser.add_option("-g", "",  dest="align_global", action="store_true",
                    help="align global")
    parser.set_defaults(align_global=False)

    parser.add_option("-s", "",  dest="align_semiglobal", action="store_true",
                    help="align semi-global")
    parser.set_defaults(align_semiglobal=False)

    parser.add_option("-p", "",  dest="gap_penalty", type="int",
                    help="Gap penalty (%default)")
    parser.set_defaults(gap_penalty=2)

        ## in the prompt line, if we write "--orf yes" then we will have the option to
## see all the possible ORFs in our sequence
    parser = optparse.OptionParser()
    parser.add_option('--vep', dest='vep',  metavar="<file>",
                        help="input alignment file (fasta)")
    parser.add_option('--vep', action='store', dest='vep', type='file')
    (options, args) = parser.parse_args()

    if options.vep == 'yes':
        vcf = vep()

    # get the options:
    (options, args) = parser.parse_args()

    if not options.fasta:
        # check if we have an option left (to be used as input filename):
        if args:
            options.fasta = args.pop()
        else:
            print ("Need at least an input file (fasta)")
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
    