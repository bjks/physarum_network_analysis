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


SHOW = False
SAVE = False
plt.rcParams["font.family"] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Avenir',
                                    'Tahoma',
                                    'DejaVu Sans',
                                    'Lucida Grande',
                                    'Verdana']

def upd_out(d, *args, supr=False):
    d.update(args)
    if not supr:
        print(*args)


def process_phase(set, label):
    align_keyword   = 'reference_point'

    times           = set.times
    positions       = set.positions

    max_period      = 140
    min_frequency   = 1./max_period

    bandwidth      = np.array([-0.002, 0.002])
    cmap_r          = 'viridis'

    data_file = branch_datfile(set, label, ext='.npz')

    if not os.path.isfile(data_file):
        print("\nNOTE: ", data_file, " not found!\n")
    else:
        print("Analyze: ", data_file)
        kymos_data = np.load(data_file)
        path_name = branch_plotpath(set, label)

        filename = path_name + '/branch_' + str(label) + '_'

        to_save_dict = dict()
        ##################################################################
        ########################## Collect data ##########################
        ##################################################################
        print("Collect data...")

        path0 = np.array(kymos_data['path'][0])

        radii = kymograph.get_dat(kymos_data, 'kymo_local_radii',
                                    'radius',
                                    set.frame_int,
                                    set.pixel_scaling,
                                    set.pixel_scaling,
                                    cmap_r,
                                    align_keyword,
                                    times=times, positions=positions)


        kymo_i = kymograph.get_dat(kymos_data, 'kymo_concentration',
                                    'intensity',
                                     set.frame_int,
                                     set.pixel_scaling,
                                     set.pixel_scaling,
                                     'Greys_r',
                                     align_keyword,
                                     times=times, positions=positions)



        print("Plotting raw...")

        plot_kymograph([radii, kymo_i], filename)


        ##################################################################
        #################### Crop Kymos #########################
        ##################################################################
        radii = radii.crop_kymo()
        kymo_i = kymo_i.crop_kymo()

        print("Plotting cropped...")
        plot_kymograph([radii, kymo_i], filename)


        print("Plotting time_series...")
        plot_time_series([radii, kymo_i], path0, filename, window_size=50)

        ##################################################################
        #################### Fourier spectrum  ###########################
        ##################################################################

        freq_r = power_spec(radii, filename, min_frequency, lowcut=0,
                            highcut=0.05, band=bandwidth,logscale=True)

        upd_out(to_save_dict, ('freq_r', freq_r) )


        band_f1, band_f2 = freq_r + bandwidth
        upd_out(to_save_dict, ('band_f1', band_f1), ('band_f2', band_f2) )

        dom_freq_in_t = spectro(radii, filename, logscale=True, window_size_in_s=300, overlap=2)

        # radii_adapt = radii.adaptive_bandpass(dom_freq_in_t, bandwidth)
        # plot_kymograph([radii_adapt], filename)

        print('Detrending...')
        radii, kymo_i = [k.detrend(200, tsig2=20, psig1=1, psig2=5,
            method = 'gauss') for k in [radii, kymo_i]]

        plot_kymograph_series([radii, kymo_i], filename)

        phase_radius, amp_radius, freq_radius, car_radius = radii.kymo_hilbert()


        plot_kymograph_series([radii, phase_radius, amp_radius, car_radius], filename,
                                ['minmax','minmax', 'minmax', 'minmax'])

        np.savez_compressed(branch_datfile(set, label, ext='_phase'), **to_save_dict)





def main():
    set_keyword = os.sys.argv[1].strip()
    color = os.sys.argv[2].strip()

    set             = data(set_keyword, no=data(set_keyword).first,
                            method='inter_mean', color=color)


    if len(os.sys.argv)>3:
        label =  int(os.sys.argv[3])
        process_phase(set, label)

    else:
        num_threads = multiprocessing.cpu_count()
        print("Number of detected cores: ", num_threads)

        labels          = get_seeds_positions(set, range_only=True)
        p = multiprocessing.Pool(num_threads)
        p.starmap(process_phase, zip(itertools.repeat(set), labels))

        p.close()
        p.join()


if __name__ == '__main__':
    main()
