# -*- coding: utf-8 -*-
import subprocess


def call_subprocess(cmd):
    # register process
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         shell=True)
    # call subprocess
    out, err = p.communicate()
    return out, err
