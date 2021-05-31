# coding: utf-8
import glob
import pandas as pd
import os


def parser(ensemblID, db_dir):
    '''
    Parse input and detect whether is a VCF or VEP file. Any other format
    is invalid. 

    Parameters
    ----------
    ensemblID : str        
        Ensemble protein id 
    db_dir : str
        directory where to find the database to parse

    Returns
    -------
    df 
        parsed file
    '''
    #sep : str
    # separation of the input file(e.g., " ", "\t", ...)

    # similar to grep. Faster than reading the
    # create empty list to store the rows
    f = glob.glob(os.path.join(db_dir, (ensemblID + '.*')))
    if not f:
        raise IOError()
        #exit(-1)
    else:
        try:
            df = pd.read_csv(f[0], sep="\t| ", engine='python')
        except:
            df = pd.read_csv(f[0], sep=" ", engine='python')
        return df
