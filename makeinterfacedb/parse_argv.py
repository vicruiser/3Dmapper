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

    parser.add_argument("--pdb", dest="pdb", nargs='+', metavar="<String>",
                        help="PDB file path", required=True)

    parser.add_argument('--blast_db', dest="blastdb", metavar="<String>",
                        help='proteome files path (output of makeblastdb)', required=True)

    # create default output directory
    parser.add_argument("-o", "--out", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out="./structural_db")

    parser.add_argument("-d", "--dist", metavar="<String>", dest="dist",
                        help="inter-residue distance threshold in angstroms (5 by default)")
    parser.set_defaults(dist=5)

    parser.add_argument("-t", "--int-type",
                        metavar="<String>",
                        dest="type",
                        choices=['noh','h', 'calpha', 'cbeta', 'sidechain', 'backbone', 'all' ],
                        help="interface definition. Options are: 'noh' (by default) to calculate distance considering heavy atoms only; 'h' to compute hydrogen bonds; \
                         'calpha' to compute CA-CA distance; 'cbeta' to measure distances between CB-CB (CA in the case of Glycine); \
                         'sidechain' consider atoms only from sidechains; 'backbone' to calculate distances between atoms in the backbone only; \
                          or 'all' to calcule interfaces between any possible atom")
    parser.set_defaults(type='noh')

    parser.add_argument("-e", "--evalue", metavar="<float>", dest="evalue",
                        help="e-value threshold in BLAST. Default is 1e-7.")
    parser.set_defaults(evalue=0.0000001)                  
    
    parser.add_argument("--pident", metavar="<float>", dest="pident",
                        help="percent identity threshold between query (pdb chain) and hit (protein) sequences. Default is 20 percent")
    parser.set_defaults(pident=20)

    parser.add_argument("-c", "--coverage", metavar="<float>", dest="coverage",
                        help="percent coverage threshold of the protein sequence (how much of the protein sequence is covered by the PDB sequence). Default is 0 percent")
    parser.set_defaults(coverage=0)

    parser.add_argument("--interaction",
                        metavar="<String>",
                        dest="int",
                        choices=['protein', 'ligand', 'nucleic', 'all'],
                        help="Interaction type: 'protein', 'ligand', 'nucleic', combination of the previous separated by space, or 'all' (by default)")
    parser.set_defaults(int='all')

    parser.add_argument("--biolip", dest="biolip", action='store_true',
                        help="consider BioLiP list (https://zhanggroup.org/BioLiP/ligand_list) to remove artifact ligands? Default is False", default=False)

    parser.add_argument('-p', "--parallel", dest="parallel", action='store_true',
                        default=False,
                        help="Parallelize process")

    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel")
    parser.set_defaults(njobs=1)

    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
