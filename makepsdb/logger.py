# -*- coding: utf-8 -*-
import time
import logging
import os


# automatize logging process
def get_logger(name, out_dir):
    log_format = '[' + '%(asctime)s' + ']' + \
        '%(name)8s -  %(message)s'
    logging.basicConfig(level=logging.DEBUG,
                        format=log_format,
                        filename=os.path.join(
                            out_dir, 'makepsdb.log'),
                        filemode='a')

    return logging.getLogger(name)
