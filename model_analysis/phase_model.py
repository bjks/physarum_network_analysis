import numpy as np
import matplotlib.pyplot as plt

import multiprocessing
import itertools

import scipy.io as sio

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

def load_mat_file(mat_file):
    mat_contents = sio.loadmat(mat_file)
    print(sorted(mat_contents.keys()))
    return mat_contents

def unique_tol(arr, tol=1e-3):
    return arr[~(np.triu(np.abs(arr[:,None] - arr) <= tol,1)).any(0)]


def phase_model(mat_file):

    times           = (100, None) # at least ~10 to get constant timeincrement
    positions       = (None, None)


    cmap_r          = 'viridis'
    cmap_c          = 'cividis'

    kymos_data = load_mat_file(mat_file)
    filename = mat_file[:-4] + '_'

    time_steps = np.diff(kymos_data['binTime'][times[0]:times[1],0])

    t_scaling = unique_tol(time_steps)
    print("time steps", np.shape(kymos_data['binTime']))
    if np.size(t_scaling)>1:
        print('WARNING: time_step not unique up to tolerance')
    else:
        t_scaling = t_scaling[0] * 10
    print('time_step: ', t_scaling)


    to_save_dict = dict()

    radii           = kymograph.get_dat(kymos_data, 'binRadius',
                                        'radius',
                                        t_scaling,
                                        1.0,
                                        100,
                                        cmap=cmap_r,
                                        align_keyword=None,
                                        unit=r'$\mu$m',
                                        times=times, positions=positions,
                                        ismat_file=True)


    mass     = kymograph.get_dat(kymos_data, 'binMass',
                                        'mass',
                                        t_scaling,
                                        1.0,
                                        1.0,
                                        cmap=cmap_c,
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True)


    flow = kymograph.get_dat(kymos_data, 'binVelocity',
                                        'flow velocity',
                                        t_scaling,
                                        1.0,
                                        1.0,
                                        cmap='Spectral',
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True)


    conce = mass.copy_meta(mass.kymo / radii.kymo**2)
    conce.name = 'concentration'

    radii.color = 'orange'
    conce.color = 'blue'
    flow.color = 'green'


    print("Plotting raw...")

    plot_kymograph([radii, mass, conce, flow], filename)


    print("Plotting time_series...")
    plot_time_series([radii, mass, conce, flow], None, filename,
                    window_size=10, spacial_avg=False)


    print('Detrending...')
    radii, mass, conce = \
        [k.detrend(cval=np.mean(k.kymo), method = 'const') \
        for k in [radii, mass, conce]]


    ##################################################################
    #################### Fourier spectrum  ###########################
    ##################################################################

    freq_r = power_spec(radii, filename, 0, lowcut=0,
                        highcut=1, band=None,logscale=True)
    freq_c = power_spec(conce, filename, 0, lowcut=0,
                        highcut=1, band=None,logscale=True)

    upd_out(to_save_dict, ('freq_r', freq_r), ('freq_c', freq_c) )


    print("Time shift map...")
    time_shift_map = time_shift2d(radii, conce, filename, upsample_t=1,
                        window_size_in_s=2./freq_r,
                        t_sampling_in_s=2/freq_r,
                        x_sampling=25, search_range_in_s=1./freq_r,
                        detrend=None)
    upd_out(to_save_dict, ('time_shift_map', time_shift_map), supr=True)



    ##################################################################
    ################# calc phase of oscillation ######################
    ##################################################################
    print('Phase...')
    phase_radius, amp_radius, freq_radius, car_radius = radii.kymo_hilbert()
    phase_conce,  amp_conce, freq_conce, car_conce = conce.kymo_hilbert()

    plot_kymograph_series([radii, phase_radius, amp_radius, car_radius],
                            filename,
                            ['minmax','minmax', 'minmax', 'minmax'])
    plot_kymograph_series([conce, phase_conce, amp_conce, car_conce],
                            filename,
                            ['minmax', 'minmax', 'minmax', 'minmax'])

    phase_shift_map = phase_shift2d(phase_radius, phase_conce, filename,
                                    t_window_in_s=1./freq_r)

    ##################################################################
    ############ bandpass kymograph before binning ##################
    ######################## Phase dependence ########################
    ##################################################################

    phase_shift_c = phase_average(phase_radius.kymo,
                                [radii, conce],
                                filename, 'concentration', no_bins=20,
                                shift_meth='fit')

    upd_out(to_save_dict, ('phase_shift_c', phase_shift_c) )

    ##################################################################
    ######################### correlations ###########################
    ##################################################################
    corr_shift_c = correlation2d(radii, conce, filename, 'kymos',
                            upsample_t=2, upsample_x=1,
                            search_range_in_s=1./freq_r)
    upd_out(to_save_dict, ('corr_shift_c', corr_shift_c) )


    corr_shift_f = correlation2d(conce, flow, filename,
                        'flow-concentration',
                        upsample_t=2, upsample_x=1,
                        search_range_in_s=1./freq_r)
    upd_out(to_save_dict, ('corr_shift_f', corr_shift_f) )


    ##################################################################
    ######################### Save in txt ############################
    ##################################################################

    np.savez_compressed(filename + 'phase', **to_save_dict)



def main():
    phase_model(os.sys.argv[1].strip())


if __name__ == '__main__':
    main()
