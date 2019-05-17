from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from correlation import *
from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter

from multiprocessing.dummy import Pool as ThreadPool

def wall_height(r, d):
    R1 = 1. - d
    hw = np.where(r < R1, np.sqrt(1 - r**2) - np.sqrt(R1**2 - r**2), np.sqrt(1 - r**2))
    return hw

def tube_height(r):
    return np.sqrt(1 - r**2)

def model_tube(r, d, enh):
    int = tube_height(r) + enh * wall_height(r, 0.1)
    return int/int[0]

def tube_background(r, I_back=0.3):
    int = tube_height(r) + I_back * (1- tube_height(r))
    return int/int[0]

def main():

    r = np.linspace(0, 1, 1000)[:-1]
    plt.plot(r, tube_height(r), c='red', label='red channel')
    plt.plot(r, wall_height(r, 0.1), c='grey', label='wall height')
    plt.plot(r, model_tube(r, 0.1, 0.5), c='green', label='green channel')
    ratio = model_tube(r, 0.1, 0.5)/tube_height(r)
    plt.plot(r, ratio/ratio[-1], c='blue', label='norm. ratio')

    plt.legend()

    plt.savefig('../model_tube.pdf')
    plt.close()


    plt.plot(r, tube_height(r), c='red', label='red channel')
    plt.plot(r, tube_background(r), c='green', label='green channel')
    plt.plot([0, r[-1]], [0.3, 0.3], c='grey', label='background intensity')
    ratio = tube_background(r)/tube_height(r)
    plt.plot(r, ratio/ratio[-1], c='blue', label='norm. ratio')

    plt.legend()
    plt.savefig('../model_background.pdf')




    plt.show()

if __name__ == '__main__':
    main()
