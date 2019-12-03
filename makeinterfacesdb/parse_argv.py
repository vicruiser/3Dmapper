# -*- coding: utf-8 -*-
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

    # filter results by type of variant
    parser.add_argument("-fv", nargs='+', dest="filter_var", metavar="<String>",
                        help="filter by variant consequence", default=None)
    parser.set_defaults(filter_var=None)

    # filter results by isoforms
    parser.add_argument("-fi", nargs='+', dest="filter_iso", metavar="<String>",
                        help="filter by APPRIS isoform", default='principal1')
    parser.set_defaults(filter_iso='principal1')

    # interfaces database file
    parser.add_argument("-force", dest="force", metavar="<String>",
                        help="force to owerwrite? (y/n)", default="y")

    # interfaces database file
    parser.add_argument("-intdb", dest="intdb", metavar="<file>",
                        help="interfaces database directory")
    parser.set_defaults(intdb=None)

    # interfaces database file
    parser.add_argument("-pident", dest="pident", metavar="<int>",
                        help="threshold of sequence identity (percertage)")
    parser.set_defaults(pident=50)

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
