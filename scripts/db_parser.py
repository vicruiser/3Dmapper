#!/usr/bin/env python3

import glob
import pandas as pd

def parser(ensemblID, db_dir, sep):
    # similar to grep. Faster than reading the
    # create empty list to store the rows
    f = glob.glob(db_dir + '/' + ensemblID + '*')[0]
    df = pd.read_csv(f, sep = sep)
    
    return df