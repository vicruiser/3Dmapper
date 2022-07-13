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
            if '.gz' in f: 
               df = pd.read_csv(f[0], sep="\t", engine='c', compression='gzip',
                   error_bad_lines=False)
            else: 
                df = pd.read_csv(f[0], sep="\t", engine='c')
        except:
            if '.gz' in f:
                df = pd.read_csv(f[0], sep=" ", engine='c',compression='gzip',
                   error_bad_lines=False)
            else: 
                df = pd.read_csv(f[0], sep=" ", engine='c')
        return df
