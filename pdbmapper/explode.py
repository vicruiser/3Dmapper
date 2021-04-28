# -*- coding: utf-8 -*-
import sys
import os
import re
import pandas as pd
import numpy as np

def explode(df,column_selectors, row_delimiter):
    # we need to keep track of the ordering of the columns
    def _split_list_to_rows(row,row_accumulator,column_selector,row_delimiter):
        split_rows = {}
        max_split = 0
        for column_selector in column_selectors:

            split_row = row[column_selector].split(row_delimiter)
            split_rows[column_selector] = split_row
            if len(split_row) > max_split:
                max_split = len(split_row)
            
        for i in range(max_split):
            new_row = row.to_dict()
            for column_selector in column_selectors:
                try:
                    new_row[column_selector] = split_rows[column_selector].pop(0)
                except IndexError:
                    new_row[column_selector] = ''
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(_split_list_to_rows,axis=1,args = (new_rows,column_selectors,row_delimiter))
    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df


# def explode(df, columns):
#      print(df)
#      print(columns)
#      idx = np.repeat(df.index, df[columns[0]].str.len())
#      a = df.T.reindex(columns).values
#      #print(a[1])
#      #print(np.concatenate(a[1]))
#      concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
#      print(concat)
#      p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
#      return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

# credit to @piRSquared from stackoverflow ###################################
# def explode(df, columns):
#     idx = np.repeat(df.index, df[columns[0]].str.len())
#     #print(idx)
#     a = df.T.reindex(columns).values
#     print(df.T)
#     concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
#     print(concat)
#     p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
#     print(p)
#     res = pd.concat([df.drop(columns, axis=1), p],
#                     axis=1).reset_index(drop=True)
#     return res
##############################################################################
# def explode(df, cols, split_on='-'):
#     """
#     Explode dataframe on the given column, split on given delimeter
#     """
#     print(cols)
#     cols_sep = list(set(df.columns) - set(cols))
#     print(cols_sep)
#     df_cols = df[cols_sep]
#     explode_len = df[cols[0]].str.split(split_on).map(len)
#     print(explode_len)
#     repeat_list = []
#     for r, e in zip(df_cols.as_matrix(), explode_len):
#         repeat_list.extend([list(r)]*e)
    
#     df_repeat = pd.DataFrame(repeat_list, columns=cols_sep)
#     print(df_repeat)
#     df_explode = pd.concat([df[col].str.split(split_on, expand=True).stack().str.strip().reset_index(drop=True)
#                             for col in cols], axis=1)
#     df_explode.columns = cols
#     print(df_explode)
#     return pd.concat((df_repeat, df_explode), axis=1)