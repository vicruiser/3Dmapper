# -*- coding: utf-8 -*-
# import necessary modules
import os


def isfile(infile):
    # check if infile is a txt file containing a list of files
    try:
        with open(infile) as list_files:
            # read just first line
            lines = list_files.readline()
            normpath = map(os.path.normpath, lines)
            abspath = map(os.path.abspath, normpath)
            if any(list(map(os.path.isfile, abspath))) is True:
                return 'list_files'
            else:
                raise IOError
    # check if infile is a txt file
    except IOError:
        normpath = os.path.normpath(infile)
        abspath = os.path.abspath(normpath)
        if os.path.isfile(abspath) is True:
            return 'is_file'
        else:
            return 'file_not_recognized'
