#!/usr/bin/python3
import argparse, os



# def dir_path(dirName):
#     '''Parse input interfaces database to put it in the right format.
#     Parameters
#     ----------
#     protID : str
#         Ensemble protein id 
#     interfacesDB_filepath : str
#         DESCRIPTION MISSING!!
#     Returns
#     -------
#     subset_interfaces_db
#         DESCRIPTION MISSING!!
#     '''
#     # Create target Directory if don't exist
#     if not os.path.exists(dirName):
#         os.mkdir(dirName)
#         print("Directory", dirName,  "created.")
#     else:    
#         print(dirName,
#               "is an existing directory. Results will be written in there.")


def parse_commandline():
    '''Parse input interfaces database to put it in the right format.
    Parameters
    ----------
    protID : str
        Ensemble protein id 
    interfacesDB_filepath : str
        DESCRIPTION MISSING!!
    Returns
    -------
    subset_interfaces_db
        DESCRIPTION MISSING!!
    '''
    description = '''

 ----------------------------------------------- Welcome to ------------------------------------------------------- 

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
        
        
------------------------  Map annotated genomic variants to protein interfaces data in 3D. ------------------------

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
                        help="Force to owerwrite? (y/n)", default = "y")
    parser.set_defaults(intdb=None)

    # interfaces database file
    parser.add_argument("-intdb", dest="intdb", metavar="<file>",
                        help="Interfaces database")
    parser.set_defaults(intdb=None)

    # protein id string
    parser.add_argument("-protid", nargs='+', metavar="<String>", dest="protid",
                        help="Ensembl protein id", required = True)

    # create chimera script to visualize the region of interest
    parser.add_argument("-chimera", action="store_true", dest="chimera",
                        help="generates chimeraX script")
    
    # create default output directory
    parser.add_argument("-out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out="./out/")
    
    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
