#!/usr/bin/python3
from halo import Halo


def tags(text_start, text_succeed, text_fail, emoji, **kwargs):

    verbose = kwargs.get('args.verbose', True)

    def my_decorator(func):

        def wrapper(*args):
            if verbose is True:
                spinner = Halo(text='Loading', spinner='dots12', color="cyan")
                spinner.start(text=text_start)
            try:
                f = func(*args)
                if verbose is True:
                    spinner.stop()
                    spinner.stop_and_persist(symbol=emoji, text=text_succeed)
                return f
            except IOError:
                if verbose is True:
                    spinner.fail(text=text_fail)
                exit(-1)
        return wrapper

    return my_decorator


def only_fail_tag(text_fail):

    def my_decorator(func):

        def wrapper(*args):
            spinner = Halo(text='Loading', spinner='dots')
            # spinner.start(text=text_start)
            try:
                func(*args)
                #spinner.stop_and_persist(symbol = emoji, text=text_succeed)
                return func(*args)
            except IOError:
                spinner.fail(text=text_fail)
                exit(-1)
        return wrapper

    return my_decorator
