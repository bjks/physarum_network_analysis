import matplotlib.pyplot as plt
from scipy import stats
import os

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *



def snr(data_sets, color):

    masks    = [np.load(set.file_dat + '.npz')['mask'].astype(bool) for set in data_sets]
    if color=='texas':
        raw = [read_file(set.file_raw2) for set in data_sets]
    elif color=='green':
        raw = [read_file(set.file_raw1) for set in data_sets]

    snr=[]
    for r, m in zip(raw, masks):
        signal  = np.true_divide(np.sum(r*m), np.sum(m))
        back    = np.true_divide(np.sum(r*np.invert(m)), np.sum(np.invert(m)))
        snr.append(signal/back)
    return snr



def main():
    set_keyword     = os.sys.argv[1].strip()
    color           = os.sys.argv[2].strip()
    step            = 10
    method          = 'inter_mean'

    set_inds = np.arange(data(set_keyword).first, data(set_keyword).last, step)
    # set_inds = np.arange(data(set_keyword).first, 10, step)


    data_sets = [data(set_keyword, i, color=color, method=method) for i in set_inds]


    dt = data(set_keyword).frame_int * step
    time = [i * dt for i in range(len(set_inds))]

    snr_g = snr(data_sets, 'green')
    snr_t = snr(data_sets, 'texas')


    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(time, snr_g, color = 'green')
    ax2.plot(time, snr_t, color = 'red')

    ax1.set_ylabel('SNR green')
    ax2.set_ylabel('SNR red')

    ax1.set_xlabel('time (s)')
    plt.savefig(data_sets[0].file_plot_set + '_snr.pdf', dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    main()
