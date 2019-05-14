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
from fft_analysis import *

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


def main():

    # r = np.linspace(0, 1, 1000)
    # plt.plot(r, wall_height(1, r, 0.1))
    # plt.plot(r, wall_model_profile(1, r, 0.1))
    # plt.plot(r, np.sqrt(1.**2 - r**2))
    # plt.plot(r, h)
    # R = np.linspace(20,100, 100)
    # d = 10.
    # plt.plot(R, [disk_mean_num(Ri, d) for Ri in R]  )
    #
    # plt.show()

    # set_keyword = os.sys.argv[1]
    # method= 'disk_mean'
    #
    # data_sets = [data(set_keyword, i, method) for i in range(data(set_keyword).first, data(set_keyword).last)]

    # set1 = 2
    # set2 = 4
    # set15= 3
    #
    # texas1 = read_file(data_sets[set1].file_raw2)
    # texas2 = read_file(data_sets[set2].file_raw2)
    # mask1 = create_mask(texas1, data_sets[set1].sigma, data_sets[set1].threshold, data_sets[set1].halo_sig)
    # mask2 = create_mask(texas2, data_sets[set2].sigma, data_sets[set2].threshold, data_sets[set2].halo_sig)
    #
    # texas_inter =  interpl_dye(texas1, texas2)*interpl_masks(mask1, mask2)
    #
    #
    # texas15 =read_file(data_sets[set15].file_raw2)
    # mask15 = create_mask(texas15, data_sets[set15].sigma, data_sets[set15].threshold, data_sets[set15].halo_sig)
    #
    # fig, axes = plt.subplots(2,2, figsize=(6, 6), sharex=True, sharey=True)
    # ax = axes.ravel()
    #
    # ax[0].imshow(texas1*mask1)
    # ax[0].set_title('t=0')
    #
    # ax[1].imshow(texas_inter)
    # ax[1].set_title('t=1 interpolated')
    #
    # ax[2].imshow(texas2*mask2)
    # ax[2].set_title('t=2')
    #
    # ax[3].imshow(texas15*mask15)
    # ax[3].set_title('t=1')
    #
    # plt.show()

    #
    #
    n= 100
    a = [np.arange(0,n) for i in range(n)]
    a = np.transpose(np.array(a))
    x = -np.cos(a+0.5)
    y = np.cos(a)
    plt.imshow(x)
    plt.show()
    fft_kymo(x)

    # phase_corr(x, y, 'x', 'y', 2, 'y_before_x', 'start', None, detrending='gauss', upsample=10, show=True)



    # ###############################
    # pool = ThreadPool(len(seed_positions))
    # pool.starmap(process_skeleton, zip(itertools.repeat(data_sets), seed_positions))
    # pool.close()
    # pool.join()
    ###############################

if __name__ == '__main__':
    main()
