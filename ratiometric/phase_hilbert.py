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


def upd_out(d, *args, supr=False):
    d.update(args)
    if not supr:
        print(*args)


def process_phase(set, label):
    align_keyword   = 'reference_point'

    times           = set.times
    positions       = set.positions

    substract_ecto  = False

    max_period      = 140
    min_frequency   = 1./max_period

    bandwidth      = np.array([-0.005, 0.005])
    cmap_r          = 'viridis'
    cmap_c          = 'cividis'

    data_file = branch_datfile(set, label, ext='.npz')

    if not os.path.isfile(data_file):
        print("\nNOTE: ", data_file, " not found!\n")
    else:
        print("\nAnalyze: ", data_file)
        kymos_data = np.load(data_file)
        path_name = branch_plotpath(set, label)

        filename = path_name + '/branch_' + str(label) + '_'

        kymo_file = mk_mising_dir(path_name + '/kymos/') + '/branch_' + str(label) + '_'
        shift_file = mk_mising_dir(path_name + '/shifts/') + '/branch_' + str(label) + '_'
        fourier_file = mk_mising_dir(path_name + '/fourier/') + '/branch_' + str(label) + '_'


        if os.path.isfile(set.file_raw1):
            plot_branch(kymos_data['branch_map'], label, path_name + 'bf_over',
                        set.pixel_scaling, set.file_raw1)
        else:
            print("\nNOTE: branch overlay plot skipped, file is missing: ", set.file_raw1)

        to_save_dict = dict()
        ##################################################################
        ########################## Collect data ##########################
        ##################################################################
        print("Collect data...")

        path0 = np.array(kymos_data['path'][0])


        radii           = kymograph.get_dat(kymos_data, 'kymo_local_radii',
                                            'radius',
                                            set.frame_int,
                                            set.pixel_scaling,
                                            set.pixel_scaling,
                                            cmap=cmap_r,
                                            unit=r'$\mu$m',
                                            align_keyword=align_keyword,
                                            times=times, positions=positions)


        kymo_c_green     = kymograph.get_dat(kymos_data, 'kymo_c_green',
                                            'green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Greens',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)

        kymo_inner_green = kymograph.get_dat(kymos_data, 'kymo_c_inner_green',
                                            'inner green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Greens',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)

        kymo_outer_green = kymograph.get_dat(kymos_data, 'kymo_c_outer_green',
                                            'outer green',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Greens',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)

        kymo_c_texas     = kymograph.get_dat(kymos_data, 'kymo_c_texas',
                                            'red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Reds',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)

        kymo_inner_texas = kymograph.get_dat(kymos_data, 'kymo_c_inner_texas',
                                            'inner red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Reds',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)

        kymo_outer_texas = kymograph.get_dat(kymos_data, 'kymo_c_outer_texas',
                                            'outer red',
                                             set.frame_int,
                                             set.pixel_scaling,
                                             1.,
                                             cmap='Reds',
                                             align_keyword=align_keyword,
                                             times=times, positions=positions)


        if set.analyse_flow:
            flow_x = kymograph.get_dat(kymos_data, 'kymo_flow_x',
                                        'flow x',
                                        set.frame_int,
                                        set.pixel_scaling,
                                        1.,
                                        cmap='Spectral',
                                        align_keyword=align_keyword,
                                        times=times, positions=positions)

            flow_y =   kymograph.get_dat(kymos_data, 'kymo_flow_y',
                                        'flow y',
                                        set.frame_int,
                                        set.pixel_scaling,
                                        1.,
                                        cmap='Spectral',
                                        align_keyword=align_keyword,
                                        times=times, positions=positions)


        print("Plotting raw...")

        plot_kymograph([radii, kymo_c_green, kymo_c_texas], kymo_file)

        if set.analyse_flow:
            plot_kymograph([flow_x, flow_y], kymo_file)


        # ==================================================================== #
        # Norm. concentration
        # ==================================================================== #
        radii = radii.shift(set.symm_setup, 'orange')

        conce = kymograph.ca_concentration(kymo_c_green, kymo_c_texas,
                                            set.symm_setup,
                                            'concentration', cmap_c, 'blue')

        inner = kymograph.ca_concentration(kymo_inner_green, kymo_inner_texas,
                                            set.symm_setup,
                                            'inner concentration', cmap_c, 'lightskyblue')

        outer = kymograph.ca_concentration(kymo_outer_green, kymo_outer_texas,
                                            set.symm_setup,
                                            'outer concentration', cmap_c, 'darkblue')



        if set.analyse_flow:
            flow_x = flow_x.crop_kymo()
            flow_x.color = 'lightgreen'

            flow_y = flow_y.crop_kymo()
            flow_y.color = 'darkgreen'


        print("Plotting cropped...")
        plot_kymograph([radii, conce, inner, outer], kymo_file)


        print("Plotting time_series...")
        plot_time_series([radii, conce, inner, outer], path0, kymo_file,
                        window_size=50)

        ##################################################################
        #################### Fourier spectrum  ###########################
        ##################################################################

        freq_r = power_spec(radii, fourier_file, min_frequency, lowcut=0,
                            highcut=0.05, band=bandwidth,logscale=True)
        freq_c = power_spec(conce, fourier_file, min_frequency, lowcut=0,
                            highcut=0.05, band=bandwidth,logscale=True)
        upd_out(to_save_dict, ('freq_r', freq_r), ('freq_c', freq_c) )

        ##################################################################
        #################### Correlation Analysis  #######################
        ##################################################################

        print('Detrending...')
        radii_d, conce_d, inner_d, outer_d = \
            [k.detrend(tsig1=200, tsig2=10, psig1=1, psig2=5,
            method = 'gauss') for k in [radii, conce, inner, outer]]
        plot_kymograph_series([radii_d, conce_d], kymo_file)
        plot_kymograph([radii_d, conce_d, inner_d, outer_d], kymo_file)


        print('Correlation1d...')
        corr_shifts_c = correlation1d(radii_d, [conce_d, inner_d, outer_d ],
                                        shift_file, 'conce',
                                        upsample_t=4,
                                        search_range_in_s=1./freq_r)

        print('Correlation2d...')
        correlation2d(radii_d, conce_d,
                        shift_file, 'conce',
                        upsample_x=1,
                        upsample_t=4,
                        search_range_in_s=1./freq_r)

        upd_out(to_save_dict, ('corr_shifts_c', corr_shifts_c) )


        if set.analyse_flow:
            flow_y_d = flow_y.detrend(tsig1=200, tsig2=10, psig1=1, psig2=5,
                                        method = 'gauss')
            corr_shift_f = correlation2d(conce_d, flow_y_d,
                                        shift_file, 'flow-concentartion',
                                        upsample_t=2,
                                        upsample_x=1,
                                        search_range_in_s=1./freq_r)

            upd_out(to_save_dict, ('corr_shift_f', corr_shift_f) )


        print("Time shift map...")
        time_shift_map = time_shift2d(radii_d, conce_d, shift_file, upsample_t=10,
                            window_size_in_s=2./freq_r,
                            t_sampling_in_s=2/freq_r,
                            x_sampling=100, search_range_in_s=1./freq_r,
                            detrend=None)

        upd_out(to_save_dict, ('time_shift_map', time_shift_map), supr=True)



        ##################################################################
        ######################## Phase Analysis ##########################
        ##################################################################

        print('Bandpass ...')
        band_f1, band_f2 = freq_r + bandwidth
        upd_out(to_save_dict, ('band_f1', band_f1), ('band_f2', band_f2) )

        dom_freq_in_t = spectro(radii, fourier_file, logscale=True,
                                window_size_in_s=400, overlap=2)

        radii_b, conce_b, inner_b, outer_b = \
            [k.bandpass(band_f1, band_f2) for k in [radii, conce, inner, outer]]

        plot_kymograph_series([radii_b, conce_b], kymo_file)
        plot_kymograph([radii_b, conce_b, inner_b, outer_b], kymo_file)

        print('Hilbert...')
        phase_radius, amp_radius, f_radius, car_radius = radii_b.kymo_hilbert()
        phase_conce,  amp_conce, f_conce, car_conce = conce_b.kymo_hilbert()

        print('phase_shift2d...')
        phase_shift_map = phase_shift2d(phase_radius, phase_conce, shift_file)

        plot_kymograph([phase_radius], kymo_file)
        plot_kymograph_series([radii, phase_radius, amp_radius, car_radius],
                                kymo_file,
                                ['minmax','minmax', 'minmax', 'minmax'])
        plot_kymograph_series([conce, phase_conce, amp_conce, car_conce],
                                kymo_file,
                                ['minmax', 'minmax', 'minmax', 'minmax'])


        print('Phase average...')
        phase_shift_c = phase_average(phase_radius.kymo,
                                    [radii_b, conce_b, inner_b, outer_b],
                                    shift_file, 'concentration (a.u.)',
                                    no_bins=20)

        upd_out(to_save_dict, ('phase_shift_c', phase_shift_c) )

        if set.analyse_flow:
            flow_x_b = flow_x.bandpass(band_f1, band_f2)
            flow_y_b = flow_y.bandpass(band_f1, band_f2)

            plot_kymograph([flow_x_b, flow_y_b], kymo_file)

            phase_shift_f = phase_average(phase_radius.kymo,
                                        [radii, flow_x_b, flow_y_b],
                                        kymo_file, 'flow', no_bins=20)

            upd_out(to_save_dict, ('phase_shift_f', phase_shift_f) )


        ##################################################################
        ######################### Save in npz ############################
        ##################################################################

        np.savez_compressed(branch_datfile(set, label, ext='_phase')
                                            , **to_save_dict)
        print('Saved npz \n')



def main():
    set_keyword     = os.sys.argv[1].strip()

    set             = data(set_keyword, no=data(set_keyword).first,
                            method='inter_mean', color='sep')

    if len(os.sys.argv)>2:
        label =  int(os.sys.argv[2])
        process_phase(set, label)

    else:
        num_threads = multiprocessing.cpu_count()
        print("Number of detected cores: ", num_threads)

        labels = get_seeds_positions(set, range_only=True)
        p = multiprocessing.Pool(num_threads)
        p.starmap(process_phase, zip(itertools.repeat(set), labels))

        p.close()
        p.join()

    log_message(set, 'phase_hilbert')

if __name__ == '__main__':
    main()
