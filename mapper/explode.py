# -*- coding: utf-8 -*-
import sys
import os
import re
import pandas as pd
import numpy as np

def explode(df,cols,split_on):
    """
    Explode dataframe on the given column, split on given delimeter
    """
    cols_sep = [x for x in list(df.columns) if x not in cols]

    df_cols = df[cols_sep]
    explode_len = df[cols[0]].str.split(split_on).map(len)

    repeat_list = []
    for r, e in zip(df_cols.values, explode_len):
        repeat_list.extend([list(r)]*e)
    
    df_repeat = pd.DataFrame(repeat_list, columns=cols_sep)

    
    df_explode = pd.concat([df[col].str.split(split_on, expand=True).stack().str.strip().reset_index(drop=True)
                            for col in cols], axis=1)
    df_explode.columns = cols
    
    return pd.concat((df_repeat, df_explode), axis=1)
