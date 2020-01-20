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

    # interfaces database file
    parser.add_argument("-force", dest="force", metavar="<String>",
                        help="force to owerwrite? (y/n)", default="y")

    # interfaces database file
    parser.add_argument("-intdb", nargs='+', metavar="<String>",
                        help="interfaces database directory", required=True)
    parser.set_defaults(intdb=None)

    # create default output directory
    parser.add_argument("-out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out=".")

    # create default output directory
    parser.add_argument("-v", dest="verbose", action='store_true',
                        default=False, help="Show progress of PDBmapper")

    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
