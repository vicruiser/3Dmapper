# -*- coding: utf-8 -*-
# import necessary modules
import os


def isfile(infile):

    try:
        normpath = os.path.normpath(infile)
        abspath = os.path.abspath(normpath)
        if os.path.isfile(abspath) is True:
            return 'yes'
        else:
            return 'no'
    except:
        return 'not_recognized'
