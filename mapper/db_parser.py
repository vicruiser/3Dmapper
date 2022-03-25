# coding: utf-8
import glob
import pandas as pd
import os
#import dask.dataframe as dd


def parser(prot_id, db_dir):
    '''
    Parse input and detect whether is a VCF or VEP file. Any other format
    is invalid. 

    Parameters
    ----------
    prot_id : str        
        protein id 
    db_dir : str
        directory where to find the database to parse

    Returns
    -------
    df 
        parsed file
    '''

    f = glob.glob(os.path.join(db_dir, (prot_id + '.*')))
    if not f:
        raise IOError()
    else:
        try:
<<<<<<< HEAD
            if '.gz' in f: 
               df = pd.read_csv(f[0], sep="\t| ", engine='python', compression='gzip',
                   error_bad_lines=False)
            else: 
                df = pd.read_csv(f[0], sep="\t| ", engine='python')
=======
            df = pd.read_csv(f[0], sep="\t| ", engine='python')
>>>>>>> d4d3a59650c6de6b30b08112a8aa0c0773363858
        except:
            if '.gz' in f:
                df = pd.read_csv(f[0], sep=" ", engine='python',compression='gzip',
                   error_bad_lines=False)
            else: 
                df = pd.read_csv(f[0], sep=" ", engine='python')
        return df
