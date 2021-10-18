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
    ----------------------------------------------------
    ____  _____                                        
   |___ \|  __ \                                       
     __) | |  | |_ __ ___   __ _ _ __  _ __   ___ _ __ 
    |__ <| |  | | '_ ` _ \ / _` | '_ \| '_ \ / _ \ '__|
    ___) | |__| | | | | | | (_| | |_) | |_) |  __/ |   
   |____/|_____/|_| |_| |_|\__,_| .__/| .__/ \___|_|   
                                | |   | |              
                                |_|   |_|              
    ----------------------------------------------------
    '''

    epilog = '''
          -----------------------------------
         |  >>>Publication link<<<<<         |
         |  victoria.ruiz.serra@gmail.com    |
          -----------------------------------
        '''
    # innit parser
    parser = argparse.ArgumentParser(epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter, description=description)
    # interfaces database file
    parser.add_argument("-psdb", "--prot_struct_db", nargs='+', metavar="<String>", dest="psdb",
                        help="protein structure database file created with makeinterfacedb", required=True)
    parser.set_defaults(psdb=None)
    
    parser.add_argument("-o", "--out", metavar="<String>", dest="out",
                        help="output directory. Default is current directory.")
    parser.set_defaults(out=".")
    
    parser.add_argument("-f", "--force", dest="force", action='store_true',
                        help="force to owerwrite? Inactive by default", default=False)
    
    parser.add_argument("-s", "--sort", dest="sort", action='store_true',
                        help="sort input file to split ", default=False)
    
    parser.add_argument("-p", "--parallel", dest="parallel", action='store_true',
                         default=False,
                         help="Speed up running time. Depends on GNU Parallel. \
                         O. Tange(2011): GNU Parallel - The Command-Line Power Tool, \
                         login: The USENIX Magazine, February 2011: 42-47.")
    
    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel", default = 1)
    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
