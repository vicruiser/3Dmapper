# -*- coding: utf-8 -*-
import time
import logging
import os

# set up the logging
# def setup_logger(logfilename):
#     logger = logging.getLogger(__name__)  # probably __name__ == "__main__"
#     logger.setLevel(logging.DEBUG)
#     formatter = logging.Formatter('%(asctime)s - %(message)s')
#     # if logfilename:
#     # write to logfile
#     handler = logging.FileHandler(logfilename)
#     # else:
#     # write to stdout
#     #    handler = logging.StreamHandler()
#     handler.setFormatter(formatter)
#     logger.addHandler(handler)
#     return logger


def get_logger(name, out_dir):
    log_format = '[' + '%(asctime)s' + ']' + \
        '%(name)8s -  %(message)s'
    logging.basicConfig(level=logging.DEBUG,
                        format=log_format,
                        filename=os.path.join(out_dir, 'pdbmapper.log'),
                        filemode='a')
    #console = logging.StreamHandler()
    # console.setLevel(logging.DEBUG)
    # console.setFormatter(logging.Formatter(log_format))
    # logging.getLogger(name).addHandler(console)
    return logging.getLogger(name)
