

# Imports


import sys
import numpy as np
import json



# JSON encoder for numpy
# From https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def log_message(msg, callback_fun = None):
    ''' Output log, either terminal or some callback taking a string as input''' 
    if callback_fun:
        callback_fun(msg)
    else:
        print(msg)


def log_progress(iteration, total, callback_fun = None):
    ''' Progress log, either terminal or some callback taking a percantage as input''' 
    
    if callback_fun:
        completed = (iteration / float(total))
        callback_fun(completed)
    else:
        
        print_progress(iteration, total)


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)

    https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
    more info: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '*' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s\r' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()
