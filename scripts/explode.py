#!/usr/bin/env python3
import sys
import os
import re
import pandas as pd
import numpy as np


# credit to @piRSquared from stackoverflow ###################################
def explode(df, columns):
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    res = pd.concat([df.drop(columns, axis=1), p],
                    axis=1).reset_index(drop=True)
    return res
##############################################################################


# def parser(input_file, ensemblID, colnames, sep):
#     # similar to grep. Faster than reading the
#     # create empty list to store the rows
#     matches = []
#     for line in input_file:
#         if ensemblID in line:
#             matches.append(line.split(sep))

#     df = pd.DataFrame(matches)
#     df.columns = colnames

#     return df

# def VEP_getter(crossref_file, geneID):
#     matches = []
#     for line in crossref_file:
#         if geneID in line:
#             vepf = line.split("\t")[1].strip()
#     return vepf
