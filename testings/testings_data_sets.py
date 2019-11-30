from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from analysis.tools import *
from analysis.kymograph_analysis import *


from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter

from multiprocessing.dummy import Pool as ThreadPool

def main():
    # set_keyword = os.sys.argv[1]
    # data_sets = [data(set_keyword, i, method='inter_mean', color='tg') for i in range(data(set_keyword).first, data(set_keyword).last)]

    k1=kymograph(1, 'radius', 1,2,3)
    k2=kymograph(2, 'conce', 1,2,3)

    kymo_list = [k1,k2]
    print(k1.name)
    print(k1.name +''.join([k.name for k in kymo_list]))
    ###############################

if __name__ == '__main__':
    main()
