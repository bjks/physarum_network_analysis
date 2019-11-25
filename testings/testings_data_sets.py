from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *



from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter
from phase_hilbert import *

from multiprocessing.dummy import Pool as ThreadPool

def main():
    set_keyword = os.sys.argv[1]
    data_sets = [data(set_keyword, i, method='inter_mean', color='tg') for i in range(data(set_keyword).first, data(set_keyword).last)]
    print(data_sets[0].times)

    ###############################

if __name__ == '__main__':
    main()
