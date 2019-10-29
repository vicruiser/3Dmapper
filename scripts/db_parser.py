#!/usr/bin/env python3
import glob
import pandas as pd
from scripts.decorator import tags

# @tags(text_start = "Parsing...",
#       text_succeed = "Parsing...done.",
#       text_fail = "Parsing...failed!",
#       emoji = "🦸")
def parser(ensemblID, db_dir, sep):
    # similar to grep. Faster than reading the
    # create empty list to store the rows
    
    f = glob.glob(db_dir + ensemblID + '*')
    
    if not f:
        raise IOError()
        exit(-1)
    else:
        df = pd.read_csv(f[0], sep = sep)
        return df
        
    