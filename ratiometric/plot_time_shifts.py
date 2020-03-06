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
                collection.append(data[key])
    return sets_labels, np.array(collection)


def plot_all_shifts(all_sets):
    plot_file = all_sets[0].file_meta_plot

    sets_labels, corr_shifts_c = collect_from_npzs(all_sets, 'corr_shifts_c', 1)

    fig, ax = plt.subplots()
    ax.scatter(sets_labels, corr_shifts_c[:,0], c='blue', label=r'$c$')
    ax.scatter(sets_labels, corr_shifts_c[:,1], c='lightskyblue', label=r'$c^{i}$')
    ax.scatter(sets_labels, corr_shifts_c[:,2], c='darkblue', label=r'$c^{o}$')

    ax.set_ylabel(r'time shift $\Delta \mathcal{T}$ (s)')
    ax.set_xlabel('data set, branch')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'time_shifts.pdf', dpi=400,
                    bbox_inches='tight')

    # plt.show()
    plt.close()


    _, freq_r = collect_from_npzs(all_sets, 'freq_r')
    fig, ax = plt.subplots()
    ax.scatter(sets_labels, corr_shifts_c[:,0]*freq_r, c='blue', label=r'$c$')
    ax.scatter(sets_labels, corr_shifts_c[:,1]*freq_r, c='lightskyblue', label=r'$c^{i}$')
    ax.scatter(sets_labels, corr_shifts_c[:,2]*freq_r, c='darkblue', label=r'$c^{o}$')

    ax.set_ylabel(r'norm. time shift $\Delta \mathcal{T}/T$')
    ax.set_xlabel('data set, branch')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'norm_time_shifts.pdf', dpi=400,
                    bbox_inches='tight')

    # plt.show()
    plt.close()



    fig, ax = plt.subplots()
    _, phase_shift_c = collect_from_npzs(all_sets, 'phase_shift_c')
    ax.scatter(sets_labels, phase_shift_c[:,1], c='blue', label=r'$c$')
    ax.scatter(sets_labels, phase_shift_c[:,2], c='lightskyblue', label=r'$c^{i}$')
    ax.scatter(sets_labels, phase_shift_c[:,3], c='darkblue', label=r'$c^{o}$')
    ax.set_ylabel(r'phase difference $\Delta \phi$ (rad)')
    ax.set_xlabel('data set, branch')
    plt.xticks(rotation=45, ha='right')

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center', ncol=3)
    plt.savefig(plot_file + 'phase_diff.pdf', dpi=400,
                bbox_inches='tight')
    # plt.show()
    plt.close()



def main():
    use_back_corr_files = int(os.sys.argv[1].strip())

    if use_back_corr_files:
        print("\n   ------------------- USE BACK CORR ------------------- \n")
        dummy_lower_thresh = 1
    else:
        print("\n   ================ DONT USE BACK CORR ================ \n")
        dummy_lower_thresh = None

    data_keys = []
    # data_keys.append('2019-07-04')
    # data_keys.append('2019-08-08')
    data_keys.append('2019-06-20-2-r')

    data_keys.append('2019-07-03-3')

    # data_keys.append('2019-08-09-3')
    data_keys.append('2019-08-11')
    data_keys.append('2019-08-25')
    # data_keys.append('2019-08-29')

    data_keys.append('2019-09-11')
    data_keys.append('2019-09-25-1')

    # data_keys.append('2019-10-01-pgol2')

    data_keys.append('2020-01-22-1')

    all_sets = [data(key, no=data(key).first, method='inter_mean',
                    color='sep') for key in data_keys]

    for set in all_sets:
        set.lower_thresh = dummy_lower_thresh
        set.update()

    plot_all_shifts(all_sets)

if __name__ == '__main__':
    main()
