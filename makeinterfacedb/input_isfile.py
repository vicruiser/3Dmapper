# -*- coding: utf-8 -*-
# import necessary modules
import os


def isfile(arg):
    # check if infile is a txt file containing a list of files
    for infile in arg: 
        try:
            with open(infile) as list_files:
                # read just first line
                lines = list_files.readline().strip('\n')
                #normpath = map(os.path.normpath, lines)
                #abspath = map(os.path.abspath, normpath)                
                if os.path.isfile(lines) is True:
                    return 'list_files'
                    break
                else:
                    raise IOError
                
        # check if infile is a txt file
        except (UnicodeDecodeError, IOError, OSError):
            normpath = os.path.normpath(infile)
            abspath = os.path.abspath(normpath)
            if os.path.isfile(abspath) is True:
                return 'is_file'
            else:
                return 'file_not_recognized'
