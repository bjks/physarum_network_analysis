import numpy as np
import matplotlib.pyplot as plt

import multiprocessing
import itertools

import os

import sys
sys.path.append("..")

from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.kymograph_analysis import *
from analysis.plotting import *


plt.rcParams["font.family"] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['CMU Sans Serif']
plt.rcParams['mathtext.default'] = 'regular'
params = {'text.usetex': False, 'mathtext.fontset': 'cm'}
plt.rcParams.update(params)


def collect_from_npzs(sets, key, print_out=False):
    collection = []
    sets_labels = []
    for set in sets:
        if print_out:
            print(set.keyword)
        for label in range(10):
            data_file = branch_datfile(set, label, ext='_phase.npz')
            if not os.path.isfile(data_file):
                pass
                # print("\nNOTE: ", data_file, " not found!\n")
            else:
                if print_out:
                    print(branch_datfile(set, label, ext='_phase.npz'))
                data = np.load(data_file)
                sets_labels.append(set.keyword + ', ' + str(label))
                collection.append(np.array(data[key]))
    return sets_labels, np.array(collection)


def plot_all_shifts(all_sets, sets_labels, full):
    plot_file = all_sets[0].file_meta_plot

    if full:
        plot_file += '_' + full + '_'

    if sets_labels == None:
        sets_labels, corr_shifts_c = collect_from_npzs(all_sets, 'corr_shifts_c', 1)
    else:
        _, corr_shifts_c = collect_from_npzs(all_sets, 'corr_shifts_c', 1)
    print(corr_shifts_c)

    size = 80
    size_small = size*0.66
    c_dir = dict(c='blue', label=r'$C$', s=size, alpha=0.5)
    ci_dir = dict(c='lightskyblue', label=r'$C^{i}$', s=size_small, alpha=0.5)
    co_dir = dict(c='darkblue', label=r'$C^{o}$', s=size_small, alpha=0.5)


    fig, ax = plt.subplots()
    ax.scatter(sets_labels, corr_shifts_c[:,0], **c_dir)
    ax.scatter(sets_labels, corr_shifts_c[:,1], **ci_dir)
    ax.scatter(sets_labels, corr_shifts_c[:,2], **co_dir)

    ax.set_ylabel(r'time shift $\Delta \mathcal{T}$ (s)')
    # ax.set_xlabel('data set')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 1.1), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'time_shifts.pdf', dpi=400,
                    bbox_inches='tight')

    # plt.show()
    plt.close()


    _, freq_r = collect_from_npzs(all_sets, 'freq_r')
    fig, ax = plt.subplots()
    ax.scatter(sets_labels, corr_shifts_c[:,0]*freq_r, **c_dir)
    ax.scatter(sets_labels, corr_shifts_c[:,1]*freq_r, **ci_dir)
    ax.scatter(sets_labels, corr_shifts_c[:,2]*freq_r, **co_dir)

    ax.set_ylabel(r'norm. time shift $f \Delta \mathcal{T}$')
    # ax.set_xlabel('data set')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 1.1), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'norm_time_shifts.pdf', dpi=400,
                    bbox_inches='tight')

    # plt.show()
    plt.close()



    fig, ax = plt.subplots()

    ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
    labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']

    _, phase_shift_c = collect_from_npzs(all_sets, 'phase_shift_c')
    print(phase_shift_c)
    ax.scatter(sets_labels, phase_shift_c[:,0], **c_dir)
    ax.scatter(sets_labels, phase_shift_c[:,1], **ci_dir)
    ax.scatter(sets_labels, phase_shift_c[:,2], **co_dir)
    ax.set_ylabel(r'phase shift $\Delta \phi$ (rad)')

    ax.set_yticklabels(labels)
    ax.set_yticks(ticks)
    # ax.set_xlabel('data set')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 1.1), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'phase_diff.pdf', dpi=400,
                bbox_inches='tight')
    # plt.show()
    plt.close()



def main():
    # use_back_corr_files = int(os.sys.argv[1].strip())
    #
    # if use_back_corr_files:
    #     print("\n   ------------------- USE BACK CORR ------------------- \n")
    #     dummy_lower_thresh = 1
    # else:
    #     print("\n   ================ DONT USE BACK CORR ================ \n")
    #     dummy_lower_thresh = None


    full = 'pp'

    data_keys = []
    sets_labels = []


    if full == 'full':
        data_keys.append('2019-08-29')
        sets_labels.append('Data set A')

        data_keys.append('2020-01-22-1')
        sets_labels.append('Data set B')

        data_keys.append('2019-06-20-2-r')
        sets_labels.append('Data set C')

        data_keys.append('2019-07-03-3')
        sets_labels.append('Data set D')

        data_keys.append('2019-08-08')
        sets_labels.append('Data set E')

        data_keys.append('2019-08-11')
        sets_labels.append('Data set F')

        data_keys.append('2019-08-25')
        sets_labels.append('Data set G')

        data_keys.append('2019-09-11')
        sets_labels.append('Data set H')

        data_keys.append('2019-09-13')
        sets_labels.append('Data set I')

        data_keys.append('2019-09-25-1')
        sets_labels.append('Data set J')


    if full == 'pp':
        data_keys.append('2019-08-29')
        sets_labels.append('Data set A')

        data_keys.append('2020-01-22-1')
        sets_labels.append('Data set B')

        data_keys.append('2019-07-03-3')
        sets_labels.append('Data set D')

        data_keys.append('2019-09-13')
        sets_labels.append('Data set I')

        data_keys.append('2019-09-25-1')
        sets_labels.append('Data set J')



    all_sets = [data(key, no=data(key).first, method='inter_mean',
                    color='sep') for key in data_keys]

    for dummy_lower_thresh in [None,1]:
        for set in all_sets:
            try:
                set.lower_thresh = dummy_lower_thresh
                set.update()
                plot_all_shifts(all_sets, sets_labels, full)

            except Exception:
                print("Missing Data for: ", dummy_lower_thresh)


if __name__ == '__main__':
    main()
