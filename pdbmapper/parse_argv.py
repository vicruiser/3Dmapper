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

    # protein id, gene id or variant id
    annovar_group = parser.add_mutually_exclusive_group(required=True)
    annovar_group.add_argument('--ensid', nargs='+', metavar="<String>", dest="ensemblid",
                               help='single or list of Ensembl protein/gene ids provided \
                                    via command line or from a file')
    annovar_group.add_argument('--varid', nargs='+', metavar="<String>", dest="varid",
                               help='single or list of variants ids provided \
                                    via command line or from a file')
    # interfaces database file
    parser.add_argument("--intdb", dest="intdb", metavar="<file>",
                        help="interfaces database directory", required=True)
    parser.set_defaults(intdb=None)

    # interfaces database file
    parser.add_argument("--uniprot", dest="uniprot", action='store_true',
                        help="protein features are defined by Uniprot IDs", default=False)

    # variants db
    parser.add_argument("--vardb", dest="vardb", metavar="<String>",
                        help='variants database directory', required=True)
    parser.set_defaults(vardb=None)

    # create default output directory
    parser.add_argument("--out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out="./pdbmapper_results")

    # filter results by type of variant
    parser.add_argument("--consequence", nargs='+', dest="consequence", metavar="<String>",
                        help="filter by consequence type, e.g.:'missense_variant'. \
                            The set of consequences is defined by Sequence Ontology (http://www.sequenceontology.org/).", default=None)
    parser.set_defaults(filter_var=None)

    # filter results by isoforms
    parser.add_argument("--isoform", nargs='+', dest="isoform", metavar="<String>",
                        help="filter by a single or a list of APPRIS isoforms. \
                             The principal isoform is set by default. \
                             Options are: principal1, principal2, ...")
    parser.set_defaults(filter_iso=None)

    # interfaces database file
    parser.add_argument("--pident", dest="pident", metavar="<int>",
                        help="threshold of sequence identity (percertage)")
    parser.set_defaults(pident=None)

    # force overwrite
    parser.add_argument("--force", dest="force", action='store_true',
                        help="force to owerwrite? (y/n)", default=False)
    # create default output directory
    parser.add_argument("--parallel", dest="parallel", action='store_true',
                        default=False,
                        help="Speed up running time. Depends on GNU Parallel. \
                        O. Tange(2011): GNU Parallel - The Command-Line Power Tool, \
                        login: The USENIX Magazine, February 2011: 42-47.")

    # interfaces database file
    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel")
    parser.set_defaults(njobs=None)

    # create chimera script to visualize the region of interest
    # parser.add_argument("-chimera", action="store_true", dest="chimera",
    #                     help="generates chimeraX script")

    # # create default output directory
    # parser.add_argument("-v", dest="verbose", action='store_true',
    #                     default=False, help="Show progress of PDBmapper")

    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
