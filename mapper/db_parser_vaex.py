# coding: utf-8
import glob
#import pandas as pd
import os
import vaex
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
            if '.gz' in f[0]:
                try:  
                    df = vaex.from_csv(f[0], sep="\t", compression='gzip', dtype = str)
                    return df
                except: 
                    df = vaex.from_csv(f[0], sep="\t", compression='gzip', low_memory=False, dtype = str)
                    return df
            else: 
                try: 
                    df = vaex.from_csv(f[0], sep="\t", dtype = str)
                    return df
                except: 
                    df = vaex.from_csv(f[0], sep="\t", low_memory=False, dtype = str)
                    return df
        except:
            if '.gz' in f[0]:
                try:
                    df = vaex.from_csv(f[0], sep=" ", compression='gzip', dtype = str)
                    return df
                except: 
                    df = vaex.from_csv(f[0], sep=" ", compression='gzip', low_memory=False, dtype = str)
                    return df
            else: 
                try: 
                    df = vaex.from_csv(f[0], sep=" ", dtype = str)
                    return df
                except: 
                    df = vaex.from_csv(f[0], sep=" ", low_memory=False, dtype = str)
                    return df
