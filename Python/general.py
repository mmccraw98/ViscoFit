from numpy import cumsum, convolve, ones, exp, diff, linspace, pi, random, tile
from time import time
import psutil
import pickle
import os
import re
from pandas import read_csv, read_excel


def deg_to_rad(angle):
    return angle * pi / 180

def rad_to_deg(angle):
    return angle * 180 / pi


def movingaverage(arr, n=10):
    ret = cumsum(arr, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def fastmovingaverage(arr, n=10):
    return convolve(arr, ones(n), 'valid') / n


def integrate(arr, dX, n=10):
    return movingaverage(arr, n=n)


def differentiate(arr, dX):
    return diff(arr) / dX


def tic():
    global current_time
    current_time = time()


def toc(return_numeric=False):
    if return_numeric:
        return time() - current_time
    else:
        print('process completed in {:.2f}s'.format(time() - current_time))


def getmemuse(return_numeric=False):
    if return_numeric:
        return psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    else:
        print('{} mb used'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2))


def logistic(x, a=1, s=0, k=1):
    return a / (1 + exp(-k * (x - s)))


def fft():
    pass


def ifft():
    pass


def regression():
    pass


def safesave(thing, path, overwrite=False):
    '''
    safely saves a thing to a given path, avoids overwrite issues
    :param thing: obj thing to be saved
    :param path: os.path path where the thing will be saved
    :param overwrite: bool determines whether or not to overwrite
    :return: none
    '''

    def getnextfile(filename):
        '''
        given a tc_data name, with an extension, find the next numeric instance of the tc_data name
        i.e. the_file1.csv -> the_file2.csv
        :param file_name: str tc_data name with a tc_data extension
        :return: str the next numeric instance of the tc_data name
        '''
        name, extension = filename.split(sep='.')  # split the tc_data name at the tc_data extension
        # \d indicates numeric digits, $ indicates the end of the string
        stripped_name = re.sub(r'\d+$', '', name)  # remove any numbers at the end of the tc_data
        fnum = re.search(r'\d+$', name)  # get any numbers at the end of the tc_data
        # if there are any numbers at the end of the tc_data, add 1 to get the next tc_data number and cast it as a string
        # if there aren't any numbers at the end of the tc_data, use 1 as the next number
        next_fnum = str(int(fnum.group()) + 1) if fnum is not None else '1'
        return stripped_name + next_fnum + '.' + extension  # return the next tc_data number

    # defining some variables that will be used often
    dir_name = os.path.dirname(path)  # directory name
    file_name = os.path.basename(path)  # tc_data name

    # check if path exists and make it if it doesn't
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)

    # if path exists and the tc_data exists get the next available tc_data name and adjust the path to reflect the change of name
    # if overwrite is enabled, then skip the renaming step and just overwrite using the given path
    while os.path.isfile(path := os.path.join(dir_name, file_name)) and not overwrite:
        file_name = getnextfile(file_name)

    # get the tc_data extension and save the tc_data accordingly
    extension = file_name.split(sep='.')[-1]
    if extension == 'csv':
        thing.to_csv(path, index=False)
    elif extension == 'xlsx':
        thing.to_xlsx(path, index=False)
    elif extension == 'txt':
        with open(path, 'w') as f:
            f.write(thing)
    else:
        with open(path, 'wb') as f:
            pickle.dump(thing, f)


def get_files(dir, req_ext=None):
    '''
    gets all the files in the given directory
    :param dir: str directory from which you want to load files from
    :param req_ext: optional str required tc_data extension
    :return: list of str names of the files in the given directory
    '''
    if req_ext is None:
        return [os.path.join(dir, f) for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
    else:
        return [os.path.join(dir, f) for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f)) and req_ext in f]


def get_folders(dir):
    '''
    gets all the folders in the given directory
    :param dir: str directory from which you want the sub-directories
    :return: list of str names of the sub-directories
    '''
    return [f.path for f in os.scandir(dir) if f.is_dir()]


def load(path, required_extension=None):
    '''
    loads data from a number of formats into python
    :param path: str path to thing being loaded in
    :param required_extension: str required extension for the file(s) to be loaded
    i.e. only load files with the required_extension
    :return: the data
    '''
    if not os.path.isfile(path):
        exit('tc_data does not exist')
    file_name = os.path.basename(path)  # tc_data name
    extension = file_name.split(sep='.')[-1]
    if extension == 'csv' or required_extension == 'csv':
        data = read_csv(path)
    elif extension == 'xlsx' or required_extension == 'xlsx':
        data = read_excel(path, engine='openpyxl')
    elif extension == 'txt' or required_extension == 'txt':
        with open(path, 'r') as f:
            data = f.read()
    elif extension == 'pkl' or required_extension == 'pkl':
        with open(path, 'rb') as f:
            data = pickle.load(f)
    else:
        exit('tc_data extension not yet supported: {}'.format(file_name))
    return data


def selectyesno(prompt):
    '''
    given a prompt with a yes / no input answer, return the boolean value of the given answer
    :param prompt: str a prompy with a yes / no answer
    :return: bool truth value of the given answer: yes -> True, no -> False
    '''
    print(prompt)  # print the user defined yes / no question prompt
    # list of understood yes inputs, and a list of understood no inputs
    yes_choices, no_choices = ['yes', 'ye', 'ya', 'y', 'yay'], ['no', 'na', 'n', 'nay']
    # use assignment expression to ask for inputs until an understood input is given
    while (choice := input('enter: (y / n) ').lower()) not in yes_choices + no_choices:
        print('input not understood: {} '.format(choice))
    # if the understood input is a no, it returns false, if it is a yes, it returns true
    return choice in yes_choices


def altspace(start, step, count, **kwargs):
    '''
    creates an evenly spaced numpy array starting at start, with a specified step size, and a given number of steps
    :param start: float starting position of array
    :param step: float step size between consecutive elements in the array
    :param count: int total number of elements in the array
    :param kwargs: any extra arguments that may be passed in the numpy array creation
    :return:
    '''
    return linspace(start, start + (step * count), count, endpoint=False, **kwargs)


def gaussian_white_noise(amplitude, shape):
    '''
    creates a gaussian white noise signal for adding experimental noise to a signal
    :param amplitude: float amplitude (2 * standard deviation) of the noise signal
    :param num_samples: tuple of ints size of the signal
    :return: (num_samples, num_signals) numpy array with the noise signal
    '''
    return random.normal(loc=0, scale=amplitude / 2, size=shape)

def row2mat(row, n): #@TODO move to helperfunctions
    '''
    stacks a row vector (numpy (m, )) n times to create a matrix (numpy (m, n)) NOTE: CAN SLOW DOWN COMPUTATION IF DONE MANY TIMES
    :param row: numpy array row vector
    :param n: int number of replications to perform
    :return: numpy matrix (m, n) replicated row vector
    '''
    # do once at the beginning of any calculation to improve performance
    return tile(row, (n, 1)).T
