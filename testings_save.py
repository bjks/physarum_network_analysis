from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter
from phase_hilbert import *

from multiprocessing.dummy import Pool as ThreadPool

def wall_height(R, r, d):
    R1 = R - d
    hw = np.where(r < R1, np.sqrt(R**2 - r**2) - np.sqrt(R1**2 - r**2), np.sqrt(R**2 - r**2))
    return hw

def disk_mean_num(R, d):
    r = np.arange(0, R, 0.01)
    mean = np.sum(wall_height(R,r,d))/R**2
    return mean


def interpl_masks(mask1, mask2):
    mask1 = mask1.astype(bool)
    mask2 = mask2.astype(bool)

    d1 = -ndi.distance_transform_edt(mask1) + ndi.distance_transform_edt(np.invert(mask1))
    d2 = -ndi.distance_transform_edt(mask2) + ndi.distance_transform_edt(np.invert(mask2))

    d = np.mean([d1, d2], axis=0)
    new_mask = np.where(d < 0, 1, 0)
    return new_mask

def interpl_dye(raw1, raw2):
    return np.mean([raw1, raw2], axis = 0)
###############


def window_view(arr):
    return 0

def main():
    params = [1, 2]
    key = ['x', 'y']
    name = {k:p for p,k in zip(params, key)}

    print(name)
    np.savez_compressed('test_save', e=2, **name)

    dat = np.load('test_save.npz')
    print(dat['e'])
    print(dat['x'])



if __name__ == '__main__':
    main()
