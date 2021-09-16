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

    # protein id, gene id or variant id
    annovar_group = parser.add_mutually_exclusive_group(required=True)
    annovar_group.add_argument('-pid', '--prot-id', nargs='+', metavar="<String>", dest="prot_id",
                               help='single or list of Ensembl protein/gene ids provided \
                                    via command line or from a file')
    annovar_group.add_argument('-vid', '--var-id', nargs='+', metavar="<String>", dest="varid",
                               help='single or list of variants ids provided \
                                    via command line or from a file')   
    # interfaces database file
    parser.add_argument("-psdb",  dest="psdb", metavar="<String>",
                        help="interfaces database directory", required=True)

    # interfaces database file
    #parser.add_argument('-u', "--uniprot", dest="uniprot", action='store_true',
    #                    help="protein features are defined by Uniprot IDs", default=False)

    # variants db
    parser.add_argument("-vdb", '--vardb', dest="vardb", metavar="<String>",
                        help='variants database directory', required=True)

    # create default output directory
    parser.add_argument("-o", "--out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out="./3dmapper_results")
    
    # force overwrite
    parser.add_argument('-ids', dest="dict_geneprot", metavar="<String>",
                        help="File that contains protein, transcripts and gene IDs.", required = True)

    # filter results by type of variant
    parser.add_argument('-c', "--consequence", nargs='+', dest="consequence", metavar="<String>",
                        help="filter by consequence type, e.g.:'missense_variant'. \
                            The set of consequences is defined by Sequence Ontology (http://www.sequenceontology.org/).", default=None)

    # filter by isoforms
    parser.add_argument('-i', "--isoform", nargs='+', dest="isoform", metavar="<String>",
                        help="filter by a single or a list of APPRIS isoforms. \
                             The principal isoform is set by default. \
                             Options are: principal1, principal2, ...")
    # filter by distance (applicable to interfaces)
    parser.add_argument('-d', "--dist", dest="dist", metavar="<float>",
                        help="threshold of maximum distance allowed in angstroms")

    # filter by sequence identity percent
    parser.add_argument("--pident", dest="pident", metavar="<int>",
                        help="threshold of sequence identity (percertage)")
    parser.set_defaults(pident=20)
    
    # filter by sequence identity percent
    parser.add_argument('-e', "--evalue", dest="evalue", metavar="<int>",
                        help="threshold of evalue", default=None)
    # force overwrite
    file_dest = parser.add_mutually_exclusive_group(required=False)
    file_dest.add_argument('-f', "--force", dest="force", action='store_true',
                        help="force to owerwrite? (y/n)", default=False)
    file_dest.add_argument('-a', "--append", dest="append", action='store_true',
                        help="Two or more calls to the program write are able to append results to the same output file.",
                        default=False)
    
    # create default output directory
    parser.add_argument('-p', "--parallel", dest="parallel", action='store_true',
                        default=False,
                        help="Parallelize process")

    # interfaces database file
    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel")
    parser.set_defaults(njobs=1)
    # interfaces database file
    parser.add_argument("-v", "--verbose", dest="verbose", action='store_true',
                        help="Print progress.", default=False)
    parser.set_defaults(njobs=None)

    # force overwrite
    parser.add_argument('-l', "--location", dest="loc", action='store_true',
                        help="Map all variants and detect their location.", default=False)
    
    # save final results in CSV format
    parser.add_argument('-csv', dest="csv", action='store_true',
                        help="Write the contained data to a CSV file.", default=False)
    
    # save final results in HDF5 format
    parser.add_argument('-hdf', dest="hdf", action='store_true',
                        help="Write the contained data to an HDF5 file using HDFStore.", default=False)


    # create chimera script to visualize the region of interest
    # parser.add_argument("-chimera", action="store_true", dest="chimera",
    #                     help="generates chimeraX script")


    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
