# -*- coding: utf-8 -*-
# import necessary modules
import os


def isfile(infile):
    # check if infile is a txt file containing a list of files
    try:
        with open(infile) as list_files:
            lines = list_files.read().splitlines()
            normpath = map(os.path.normpath, lines)
            abspath = map(os.path.abspath, normpath)
            if any(list(map(os.path.isfile, abspath))) is True:
                return 'list_files'
            else:
                raise IOError()
    # check if infile is a txt file
    except:
        normpath = os.path.normpath(infile)
        abspath = os.path.abspath(normpath)
        if os.path.isfile(abspath) is True:
            return 'is_file'
        else:
            raise IOError()


with open(ids) as f:
                        lines = f.read().splitlines()
                        # set variable input as file
                        input = "file"
                except:
                    # set variable input as not file
                    input = "not_file"
