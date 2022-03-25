# -*- coding: utf-8 -*-
import subprocess
import os
import os.path
import re
import time
import glob
from halo import Halo
from .decorator import tags
from .logger import get_logger
from .parse_argv import parse_commandline

# add decorator to main function
@tags(text_start="Split file by selected ensembl id...This might take up some time...",
      text_succeed="Split file by selected ensembl id...done.",
      text_fail="Split file by selected ensembl id...failed!",
      emoji="\U00002702")
def index(input_file, out_dir, log_dir):
    '''
    Index variants file.

    Parameters
    ----------
    prefix : str
        Ensembl ID prefix. Either ESNG or ENSP.
    input_file : str
        Path to infile.
    out_dir : str
        Path to output.
    out_extension : str
        Output filename extension

    Returns
    -------
    indexed_file
        1st column is variants ids and the 2nd the corresponding gene ids
    '''
    # log file
    logger = get_logger('create index', log_dir)
    logger.info('Creating index file.')