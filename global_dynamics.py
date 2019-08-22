import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.data_sets import *
from scipy import stats
import os


def sum_images(data_sets, keyword):
    return np.array([np.sum(np.load(set.file_dat + '.npz')[keyword]) for set in data_sets])

def avg_images(data_sets, keyword):
    return np.array([np.mean(np.load(set.file_dat + '.npz')[keyword]) for set in data_sets])



def main():
    set_keyword = os.sys.argv[1]
    step        = int(os.sys.argv[2])
    method      = 'inter_mean'

    # set_inds = np.arange(data(set_keyword).first, data(set_keyword).last, step)
    set_inds = np.arange(data(set_keyword).first, 200, step)

    data_sets = [data(set_keyword, i, color='tg', method=method) for i in set_inds]

    area            = sum_images(data_sets, 'mask')
    # volume          = sum_images(data_sets, 'texas_clean')
    # calcium         = sum_images(data_sets, 'green_clean')
    # concentration   = sum_images(data_sets, 'ratio')/area


    dt = data(set_keyword).frame_int
    time = [i * dt for i in range(len(data_sets))]


    plt.plot(time, area/area[0], label='total area')
    # plt.plot(time, volume/volume[0], label='total volume')
    # plt.plot(time, calcium/calcium[0], label='total labeled Ca')

    # plt.plot(time, concentration/concentration[0], label='mean concentration')

    plt.xlabel('time (s)')
    plt.legend()
    plt.savefig(data_sets[0].file_plot_set + '_global_dynamics.pdf')


if __name__ == '__main__':
    main()
