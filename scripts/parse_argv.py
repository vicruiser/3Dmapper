#!/usr/bin/python3
import argparse
import os


def parse_commandline():
    '''
    Parse inputs from command line.

    Returns
    -------
    args
        arguments to give to the functions
    '''
    description = '''
----------------------------------------- Welcome to ----------------------------------------------

   $$$$$$$\  $$$$$$$\  $$$$$$$\   
   $$  __$$\ $$  __$$\ $$  __$$\     
   $$ |  $$ |$$ |  $$ |$$ |  $$ |$$$$$$\$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\ 
   $$$$$$$  |$$ |  $$ |$$$$$$$\ |$$  _$$  _$$\  \____$$\ $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$\ 
   $$  ____/ $$ |  $$ |$$  __$$\ $$ / $$ / $$ | $$$$$$$ |$$ /  $$ |$$ /  $$ |$$$$$$$$ |$$ |  \__|
   $$ |      $$ |  $$ |$$ |  $$ |$$ | $$ | $$ |$$  __$$ |$$ |  $$ |$$ |  $$ |$$   ____|$$ |
   $$ |      $$$$$$$  |$$$$$$$  |$$ | $$ | $$ |\$$$$$$$ |$$$$$$$  |$$$$$$$  |\$$$$$$$\ $$ |
   \__|      \_______/ \_______/ \__| \__| \__| \_______|$$  ____/ $$  ____/  \_______|\__|
                                                         $$ |      $$ |
                                                         $$ |      $$ |
                                                         \__|      \__|

---------------  Map annotated genomic variants to protein interfaces data in 3D. -----------------

'''
    epilog = \
        '''
          -------------------------------------------------------------------------        
         |  Copyright (c) 2019 Victoria Ruiz --                                    |  
         |  vruizser@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          -------------------------------------------------------------------------

        '''
    # innit parser
    parser = argparse.ArgumentParser(epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter, description=description)

    # vcf file
    annovar_group = parser.add_mutually_exclusive_group()
    annovar_group.add_argument('-vcf', nargs='+', metavar="<String>", dest="vcf",
                               help='input directory containing single or \
                                    multiple annotated variants files (vcf)')
    annovar_group.add_argument("-vep",   metavar="<file>",  dest="vep",
                               help="default VEP input")
    annovar_group.add_argument("-varmap", action='store_true', dest="varmap",
                               help="use ClinVar db of annotated variants",
                               default="varmap")
    #parser.set_defaults(annovar = "varmap")

    # interfaces database file
    parser.add_argument("-force", dest="force", metavar="<String>",
                        help="force to owerwrite? (y/n)", default="y")
    parser.set_defaults(intdb=None)

    # interfaces database file
    parser.add_argument("-intdb", dest="intdb", metavar="<file>",
                        help="interfaces database directory")
    parser.set_defaults(intdb=None)

    # interfaces database file
    parser.add_argument("-pident", dest="pident", metavar="<int>",
                        help="threshold of sequence identity (percertage)")
    parser.set_defaults(pident=0)

    # protein id string
    parser.add_argument("-protid", nargs='+', metavar="<String>", dest="protid",
                        help="Ensembl protein id", required=True)

    # create chimera script to visualize the region of interest
    parser.add_argument("-chimera", action="store_true", dest="chimera",
                        help="generates chimeraX script")

    # create default output directory
    parser.add_argument("-out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out="./out")

    # create default output directory
    parser.add_argument("-v", dest="verbose", action='store_true',
                        default=False, help="Show progress of PDBmapper")

    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
