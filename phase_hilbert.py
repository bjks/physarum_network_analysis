import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.signal import hilbert
from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit

import os
import itertools

from scipy import signal
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks

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



def main():
    set_keyword     = os.sys.argv[1].strip()
    color           = 'sep'
    method          = 'inter_mean'

    set             = data(set_keyword, no=data(set_keyword).first,
                            method=method, color=color)


    align_keyword   = 'reference_point'

    times           = None, None
    positions       = None, None

    substract_ecto  = False

    max_period      = 140
    range_freq      = 0.001, 0.003
    cmap_r          = 'viridis'
    cmap_c          = 'cividis'



    labels          = get_seeds_positions(set, range_only=True)

    for label in labels:
        kymos_data = np.load(brach_datfile(set, label))
        path_name = branch_plotpath(set, label)

        filename = path_name + '/branch_' + str(label) + '_'


        ##################################################################
        ########################## Collect data ##########################
        ##################################################################
        print("Collect data...")

        radii           = kymograph.get_dat(kymos_data, 'kymo_local_radii',
                                            'radius',
                                            set.frame_int,
                                            set.pixel_scaling,
                                            set.pixel_scaling,
                                            cmap_r,
                                            align_keyword,
                                            times=times, positions=positions)


        kymo_c_green     = kymograph.get_dat(kymos_data, 'kymo_c_green',
                                            'green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Greens',
                                             align_keyword,
                                             times=times, positions=positions)

        kymo_inner_green = kymograph.get_dat(kymos_data, 'kymo_c_inner_green',
                                            'inner green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Greens',
                                             align_keyword,
                                             times=times, positions=positions)

        kymo_outer_green = kymograph.get_dat(kymos_data, 'kymo_c_outer_green',
                                            'outer green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Greens',
                                             align_keyword,
                                             times=times, positions=positions)

        kymo_c_texas     = kymograph.get_dat(kymos_data, 'kymo_c_texas',
                                            'red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Reds',
                                             align_keyword,
                                             times=times, positions=positions)

        kymo_inner_texas = kymograph.get_dat(kymos_data, 'kymo_c_inner_texas',
                                            'inner red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Reds',
                                             align_keyword,
                                             times=times, positions=positions)

        kymo_outer_texas = kymograph.get_dat(kymos_data, 'kymo_c_outer_texas',
                                            'outer red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             'Reds',
                                             align_keyword,
                                             times=times, positions=positions)


        if set.analyse_flow:
            flow_x = kymograph.get_dat(kymos_data, 'kymo_flow_x',
                                        'flow x',
                                        set.frame_int,
                                        set.pixel_scaling,
                                        1.,
                                        'Spectral',
                                        align_keyword,
                                        times=times, positions=positions)

            flow_y =   kymograph.get_dat(kymos_data, 'kymo_flow_y',
                                        'flow y',
                                        set.frame_int,
                                        set.pixel_scaling,
                                        1.,
                                        'Spectral',
                                        align_keyword,
                                        times=times, positions=positions)


        print("Plotting raw...")

        plot_kymograph([radii, kymo_c_green, kymo_c_texas], filename)

        if set.analyse_flow:
            plot_kymograph([flow_x, flow_y], filename)





        ##################################################################
        #################### Norm. concentration #########################
        ##################################################################
        radii = radii.shift(set.symm_setup, 'orange')


        conce = kymograph.ca_concentration(kymo_c_green, kymo_c_texas,
                                            set.symm_setup,
                                            'concentration', cmap_c, 'blue')

        inner = kymograph.ca_concentration(kymo_inner_green, kymo_inner_texas,
                                            set.symm_setup,
                                            'inner', cmap_c, 'lightskyblue')

        outer = kymograph.ca_concentration(kymo_outer_green, kymo_outer_texas,
                                            set.symm_setup,
                                            'outer', cmap_c, 'darkblue')




        if set.analyse_flow:
            flow_x = flow_x.crop_kymo()
            flow_x.color = 'lightgreen'

            flow_y = flow_y.crop_kymo()
            flow_y.color = 'darkgreen'



        print("Plotting cropped...")

        plot_kymograph([radii, conce, inner, outer], filename)

        plot_time_series([radii, conce], filename, window_size=50)

        ##################################################################
        #################### Fourier spectrum  ###########################
        ##################################################################

        freq_r = dominant_freq(radii, max_period=max_period)
        freq_c = dominant_freq(conce, max_period=max_period)
        print('Dominant frequency(radius): ', freq_r)
        print('Dominant frequency(conce): ', freq_c)

        band_f1, band_f2 = freq_r-range_freq[0], freq_r+range_freq[1]

        power_spec(radii, filename, min_freq=0, max_freq=0.05,
                     mark_freq = freq_r, logscale=True)
        power_spec(conce, filename, min_freq=0, max_freq=0.05,
                     mark_freq = freq_c, logscale=True)


        #
        # ####################################################
        # ########### substract  Ca_global(R) ################
        # ####################################################
        # if substract_ecto:
        #     conce = substract_ecto_contribution(radii, conce, set.frame_int)
        #     inner = substract_ecto_contribution(radii, inner, set.frame_int)
        #     outer = substract_ecto_contribution(radii, outer, set.frame_int)
        #

        ##################################################################
        ################# calc phase of oscillation ######################
        ##################################################################

        phase_radius, _, _ = extract_phase(radii, filename, band_f1, band_f2)
        phase_conce, _, _  = extract_phase(conce, filename, band_f1, band_f2)


        ##################################################################
        ############ bandpass kymograph before binning ##################
        ######################## Phase dependence ########################
        ##################################################################
        radii_b, conce_b, inner_b, outer_b = \
             [k.bandpass(band_f1, band_f2) for k in [radii, conce, inner, outer]]


        phase_shift_map = phase_shift2d(radii_b, conce_b, filename, upsample=2,
                            window_size=100, sampling=100, search_range_in_s=100)

        plot_kymograph([radii_b, conce_b, inner_b, outer_b], filename)


        phase_shift = phase_average(phase_radius,
                                    [radii_b, conce_b, inner_b, outer_b],
                                    filename, 'concentration')

        if set.analyse_flow:
            flow_x_b = flow_x.bandpass(band_f1, band_f2)
            flow_y_b = flow_y.bandpass(band_f1, band_f2)

            plot_kymograph([flow_x, flow_y], filename)

            phase_shift = phase_average(phase_radius, [radii_b, flow_x, flow_y],
                                        filename, 'flow')

        ##################################################################
        ######################### correlations ###########################
        ##################################################################
        min_k = correlate_phase(radii, conce, filename, 'kymos',
                                upsample_t=10, upsample_x=1,
                                detrend='gauss')


        show_im(phase_shift_map)

        if set.analyse_flow:
            _ = correlate_phase(conce, flow_y, filename, 'flow-concentartion',
                                upsample_t=10, upsample_x=1,
                                detrend='gauss')


        ##################################################################
        ######################### Save in txt ############################
        ##################################################################
        if SAVE:
            if not os.path.exists('time_shifts_data_sets'):
                os.mkdir('time_shifts_data_sets')

            data_sets_summary = 'time_shifts_data_sets/time_shift_' + color +'.txt'
            if not os.path.isfile(data_sets_summary):
                with open(data_sets_summary, "w") as out_var:
                    out_var.write('# data_set \t label \t' +
                                    'phase_radius \t phase_conc \t' +
                                    'phase_inner \t phase_outer \t' +
                                    'min_kymo + \n')

            with open(data_sets_summary, "a") as out_var:
                out_var.write(set_keyword +
                                '\t' + str(label) +
                                '\t' + str(phase_shift[0]) +
                                '\t' + str(phase_shift[1]) +
                                '\t' + str(phase_shift[2]) +
                                '\t' + str(phase_shift[3]) +
                                '\t' + str(min_k) +'\n')


if __name__ == '__main__':
    main()
