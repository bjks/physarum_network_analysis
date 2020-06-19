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

def unique_tol(arr, tol=1e-4):
    return arr[~(np.triu(np.abs(arr[:,None] - arr) <= tol,1)).any(0)]


def phase_model(mat_file):

    kymos_data = load_mat_file(mat_file)
    filename = mat_file[:-4] + '_'

    if len(os.sys.argv)>2:
        t0 = int(os.sys.argv[2].strip())
        filename = mk_mising_dir(filename + os.sys.argv[2].strip()+ '/' ) \
                    + os.sys.argv[2].strip()
    else:
        t0 = 200. # == t0 sec

    cmap_r          = 'viridis'
    cmap_c          = 'cividis'

    #avg time_stwp of last 100 steps
    t_scaling = np.mean(np.diff(kymos_data['binTime'][-100:,0]))
    print('Total no of steps', np.shape(kymos_data['binTime']))

    #correct for unit
    t_bar = kymos_data['time_scale'][0][0]
    print('time unit:', t_bar)
    t_scaling = t_scaling * t_bar
    print('time_step: ', t_scaling)

    #correct t0

    times = (int(t0/t_scaling), None)
    positions = (None, None)

    # check wether time step is unique in time window
    time_steps = np.diff(kymos_data['binTime'][times[0]:times[1],0])
    time_steps = unique_tol(time_steps)
    if np.size(t_scaling)>1:
        print('WARNING: time_step not unique up to tolerance')
        print(t_scaling)
        return


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
                                        ismat_file=True,
                                        t0=t0)

    print(t0)

    mass     = kymograph.get_dat(kymos_data, 'binMass',
                                        'mass',
                                        t_scaling,
                                        1.0,
                                        1.0,
                                        cmap=cmap_c,
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True,
                                        t0=t0)


    flow = kymograph.get_dat(kymos_data, 'binVelocity',
                                        'flow velocity',
                                        t_scaling,
                                        1.0,
                                        1e4/t_scaling,
                                        unit=r'$\mu$m/s',
                                        cmap='Spectral',
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True,
                                        t0=t0)


    elasticity = kymograph.get_dat(kymos_data, 'binElasticity',
                                        'passive stress',
                                        t_scaling,
                                        1.0,
                                        1.0,
                                        unit=r'$E$',
                                        cmap='Reds',
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True,
                                        t0=t0)

    tension = kymograph.get_dat(kymos_data, 'binTension',
                                        'active stress',
                                        t_scaling,
                                        1.0,
                                        1.0,
                                        unit=r'$E$',
                                        cmap='Reds',
                                        align_keyword=None,
                                        times=times, positions=positions,
                                        ismat_file=True,
                                        t0=t0)


    print("shape: " ,radii.kymo.shape)
    stress = elasticity.copy_meta(tension.kymo + elasticity.kymo)
    stress.name = 'total stress'
    stress.cmap = 'Reds'

    conce = mass.copy_meta(mass.kymo / radii.kymo**2 * 100**2)
    conce.name = 'concentration'

    radii.color = 'orange'
    conce.color = 'blue'
    flow.color = 'green'
    stress.color = '#ff0000'
    tension.color = '#b20000'
    elasticity.color = '#4c0000'

    print("Plotting raw...")

    plot_kymograph([radii, mass, conce, flow, elasticity, tension, stress], filename)


    print("Plotting time_series...")
    x_lenghts = np.shape(radii.kymo)[0]
    dummy_path = np.stack([np.arange(x_lenghts)[::-1], np.zeros(x_lenghts)], axis=-1)

    plot_time_series([radii, conce, flow], dummy_path, filename,
                    window_size=33, spacial_avg=False, line=True,figsize=(14,3.2*3))

    plot_time_series([radii, conce, elasticity, tension, stress],
                    dummy_path, filename,
                    window_size=33, spacial_avg=False, filename_suf='stress',
                    figsize=(14,16), line=True)


    print('Detrending...')
    radii, mass, conce, elasticity, tension, stress = \
        [k.detrend(cval=np.mean(k.kymo), method = 'const') \
        for k in [radii, mass, conce, elasticity, tension, stress]]


    ##################################################################
    #################### Fourier spectrum  ###########################
    ##################################################################

    freq_r = power_spec(radii, filename, 0, lowcut=0,
                        highcut=0.05, band=None,logscale=True)
    freq_c = power_spec(conce, filename, 0, lowcut=0,
                        highcut=0.05, band=None,logscale=True)

    upd_out(to_save_dict, ('freq_r', freq_r), ('freq_c', freq_c) )


    print("Time shift map...")
    # time_shift_map = time_shift2d(radii, conce, filename, upsample_t=1,
    #                     window_size_in_s=2./freq_r,
    #                     t_sampling_in_s=2/freq_r,
    #                     x_sampling=25, search_range_in_s=1./freq_r,
    #                     detrend=None)
    # upd_out(to_save_dict, ('time_shift_map', time_shift_map), supr=True)



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
    phase_shift_s = phase_average(phase_radius.kymo,
                                [radii, elasticity, tension],
                                filename, 'stress terms', no_bins=20,
                                shift_meth='fit')
    phase_shift_s = phase_average(phase_radius.kymo,
                                [radii, stress],
                                filename, 'total stress', no_bins=20,
                                shift_meth='fit')

    upd_out(to_save_dict, ('phase_shift_c', phase_shift_c) )


    spectro(radii, filename, logscale=True, window_size_in_s=500, overlap=2,
                 lowcut=0., highcut=0.04)
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
