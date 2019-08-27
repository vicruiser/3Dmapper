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

    #parser.add_option("-b",  dest="fasta2", metavar="<file>",
    #                 help="input alignment file (fasta)")
    parser.set_defaults(fasta2=None)

    # get the options:
    (options, args) = parser.parse_args()

    if not options.fasta:
        # check if we have an option left (to be used as input filename):
        if args:
            options.fasta = args.pop()
        else:
            print "Need at least an input file (fasta)"
            print ""
            parser.print_help()
            print ""
            print "ERROR: no input file given"
            exit(-1)

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))

    # clean up (recommended):
    del(parser)
    return options
    