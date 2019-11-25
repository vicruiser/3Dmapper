# -*- coding: utf-8 -*-
from halo import Halo
import requests
import argparse


def tags(text_start, text_succeed, text_fail, emoji):

    def my_decorator(func):

        def wrapper(*args):

            spinner = Halo(text='Loading', spinner='dots12', color="cyan")
            spinner.start(text=text_start)
            try:
                f = func(*args)
                spinner.stop()
                spinner.stop_and_persist(symbol=emoji, text=text_succeed)
                return f
            except IOError:
                spinner.fail(text=text_fail)
                exit(-1)
        return wrapper

    return my_decorator
