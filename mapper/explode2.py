# -*- coding: utf-8 -*-
import sys
import os
import re
import pandas as pd
import numpy as np
# credit to @piRSquared from stackoverflow ###################################
def explode2(df, columns):
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    res = pd.concat([df.drop(columns, axis=1), p],
                    axis=1).reset_index(drop=True)
    return res