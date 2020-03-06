import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import rc
# from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
#                                AutoMinorLocator)
# rc('text', usetex=True)

from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit
from scipy import signal
from scipy import LowLevelCallable

import numba
from numba import cfunc, carray
from numba.types import intc, CPointer, float64, intp, voidptr

import os
import itertools
import copy

from skimage.util.shape import view_as_windows
from skimage import feature
from skimage import transform
from skimage.util.shape import view_as_windows
from skimage.feature import match_template

from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *



class kymograph:
    """
    class that collects kymograph data array, axis scaling and useful
    information for plotting:

        self.kymo       data array which is scaled (multiplied) with c_scaling

        self.name       name of the quantity eg 'radius',
        self.state      denotes current state of the analysis

        self.t_scaling  scalings used for axis labels
        self.x_scaling  ...
        self.c_scaling  ...

        self.cmap       cmap used for plotting
        self.unit       unit of data (eg 'a.u.')
        self.color      color used for plotting


    transformation methods (hilbert, bandpass, etc.) are implemented as
    constructor like methods, ie return new instances of kymograph

    methods for plotting are non-class methods


    """

    # ===================================================== #
    #  Constructor (-like) and copy methods  #
    # ===================================================== #
    def __init__(self, kymo, name, t_scaling, x_scaling, c_scaling,
                state='raw', cmap='viridis', unit='(a.u.)', color=None):

        self.kymo       = kymo * c_scaling

        self.name       = name
        self.state      = state

        self.t_scaling  = t_scaling
        self.x_scaling  = x_scaling
        self.c_scaling  = c_scaling

        self.cmap       = cmap
        self.unit       = unit
        self.color      = color

    def __getitem__(self, index):
        return self.kymo[index]

    def __setitem(self, index, dat):
        self.kymo[index] = dat



    def copy(self):
        return copy.deepcopy(self)


    def copy_meta(self, new_kymo):
        new = copy.deepcopy(self)
        new.kymo = new_kymo
        return new


    def get_file_title(self):
        return self.name.replace(' ', '') + '(' +  \
                self.state.replace(' ', '')+')'


    def get_eq(self, dep=r'$(t,x)$'):
        if self.name=='radius':
            return r'$a$' + dep
        elif self.name=='concentration':
            return r'$c$' + dep
        elif self.name.startswith('inner'):
            return r'$c^i$' + dep
        elif self.name.startswith('outer'):
            return r'$c^o$' + dep
        else:
            return ''

    def get_plot_title(self):
        if self.state == 'raw':
            return r'{0} {1} ({2})'.format(self.name, self.get_eq(), self.unit)

        else:
            return r'{0} {1} {2} ({3})'.format(self.state, self.name, self.get_eq(), self.unit)


    # ===================================================== #
    # Operations on existing data, returning new instances
    # of kymograph with state 'raw', 'cropped' or 'band'
    # ===================================================== #

    def get_dat(kymos_data, keyword, name,
                t_s=1, x_s=1, c_s=1,
                cmap='viridis', unit='a.u.',
                align_keyword='reference_point',
                times=(None,None), positions=(None, None),
                ismat_file=False):


        kymo      = kymos_data[keyword]

        if ismat_file:
            kymo = kymo[0,0]

        if align_keyword == None:
            kymo = np.transpose(kymo)
        else:
            alignment = kymos_data['alignment']
            kymo = align_kymo(kymo, align_keyword, alignment=alignment)
            kymo = np.transpose(kymo)

        kymo = kymo[slice(*positions), slice(*times)]


        return kymograph(kymo, name, t_s, x_s, c_s,
                        state='raw',
                        cmap=cmap,
                        unit=unit,
                        color=None)


    def crop_kymo(self):
        crop = crop_nans(self.kymo)

        new = self.copy_meta(crop)
        new.state = 'cropped'
        return new



    def bandpass(self, lowcut=None, highcut=None, order=5):

        band = butter_bandpass2d(self.kymo, lowcut, highcut,
                                        1./self.t_scaling, order)

        new = self.copy_meta(band)
        if lowcut == None:
            new.state = 'lowpassed'
        elif highcut == None:
            new.state = 'highpassed'
        else:
            new.state = 'bandpassed'
        return new


    def adaptive_bandpass(self, dom_freq, bandwidth, min_lowcut=0.001, order=5):

        band = adaptive_butter_2d(self.kymo, 1./self.t_scaling, dom_freq,
                                    bandwidth, order, min_lowcut)

        new = self.copy_meta(band)
        new.state = 'adaptive bandpassed'
        return new


    def detrend(self, tsig1=None, tsig2=None, psig1=None, psig2=None, cval=None,
                method='gauss'):

        if method == 'window':
            window_size = int(tsig1/self.t_scaling)
            glob = slid_aver2d(self.kymo, window_size, axis=1)

        elif method == 'gauss':
            sigma = (psig1, tsig1/self.t_scaling)
            glob = ndi.gaussian_filter(self.kymo, sigma)

        elif method == 'const':
            glob = np.ones_like(self.kymo) * cval

        detrended = self.kymo - glob

        if psig2!=None and tsig2!=None:
            sigma = (psig2, tsig2 / self.t_scaling)
            detrended = ndi.gaussian_filter(detrended, sigma)

        new = self.copy_meta(detrended)
        new.state = 'detrended'
        return new


    def kymo_hilbert(self, lowcut=None, highcut=None):

        if lowcut!=None and highcut!=None:
            kymo_b = self.bandpass(lowcut, highcut)
        else:
            kymo_b = self.copy()

        analytic_signal = signal.hilbert(kymo_b.kymo)

        amplitude_envelope = np.abs(analytic_signal)
        amp = self.copy_meta(amplitude_envelope)
        amp.state = 'amplitude'
        amp.cmap = 'plasma'

        instantaneous_phase = np.angle(analytic_signal)
        phase = self.copy_meta(instantaneous_phase)
        phase.state = 'phase'
        phase.cmap = 'twilight_shifted'
        phase.unit = 'rad'

        f_scal = 2*np.pi * kymo_b.t_scaling
        instantaneous_frequency = np.diff(np.unwrap(instantaneous_phase))/f_scal
        freq = self.copy_meta(instantaneous_frequency)
        freq.state = 'freq'
        freq.cmap = 'cool'
        freq.unit = 'Hz'


        car = np.cos(instantaneous_phase)
        carrier = self.copy_meta(car)
        carrier.state = 'carrier'
        carrier.cmap = 'cool'
        carrier.unit = 'a.u.'

        return phase, amp, freq, carrier



    # ===================================================== #
    # Create kymograph that is used in signal analysis,
    # includes cropping and shifting if setup is symm(etric) #
    # ===================================================== #
    def shift(self, symm_setup, color):
        new = self.crop_kymo()

        if symm_setup:
            shifted = ndi.zoom(new.kymo, (1,2), order=5)[1:]
            new = new.copy_meta(shifted)
            new.t_scaling = new.t_scaling/2.


        new.color = color
        return new


    def ca_concentration(green, texas, symm_setup, name, cmap, color):
        green_c = green.crop_kymo()
        texas_c = texas.crop_kymo()

        if symm_setup:
            g2 = ndi.zoom(green_c.kymo, (1,2), order=5)[:-1]
            t2 = ndi.zoom(texas_c.kymo, (1,2), order=5)[1:]

            new = green_c.copy_meta(g2/t2)
            new.t_scaling /= 2.

        else:
            new = green_c.copy_meta(green_c.kymo/texas_c.kymo)

        new.name = name
        new.cmap = cmap
        new.color = color
        return new



# ===================================================== #
#  Plotting  #
# ===================================================== #

def plot_kymograph_series(kymo_list, filename, lim_modes='min'):
    fig, axes = plt.subplots(len(kymo_list), 1)
    ax = axes.ravel()

    if np.size(lim_modes)==1:
        lim_modes = [lim_modes for i in kymo_list]

    for kymo, ax, lim_mode in zip(kymo_list, axes, lim_modes):
        extent =    [0, kymo.kymo.shape[1]*kymo.t_scaling,
                     0, kymo.kymo.shape[0]]


        if lim_mode.startswith('std'):
            stds = float(lim_mode[3:])
            std_im, mean_im = np.nanstd(kymo.kymo), np.nanmean(kymo.kymo)
            min, max = mean_im - stds*std_im, mean_im + stds*std_im

            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, extent=extent,
                            origin='lower', aspect='auto',
                            vmin = min, vmax = max)
        else:
            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, extent=extent,
                            origin='lower', aspect='auto')


        if not ax==axes[-1]:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('time (s)')


        ax.set_title(kymo.get_plot_title())
        ax.set_ylabel('space (pixel)')


        divider = make_axes_locatable(ax)
        fig.colorbar(im, cax = divider.append_axes('right', size='5%',
                    pad=0.05), orientation='vertical')

    plt.subplots_adjust(hspace=0.7)

    plt.savefig(filename + kymo_list[0].state + \
                ''.join([k.name for k in kymo_list]) + '.pdf',
                dpi=400, bbox_inches='tight')
    plt.close()



def plot_kymograph(kymo_list, filename, lim_modes='minmax'):

    if np.size(lim_modes)==1:
        lim_modes = [lim_modes for i in kymo_list]

    for kymo, lim_mode in zip(kymo_list, lim_modes):

        fig, ax = plt.subplots()

        extent =    [0, kymo.kymo.shape[1]*kymo.t_scaling,
                     0, kymo.kymo.shape[0]]

        if lim_mode.startswith('std'):
            stds = float(lim_mode[3:])
            std_im, mean_im = np.nanstd(kymo.kymo), np.nanmean(kymo.kymo)
            min, max = mean_im - stds*std_im, mean_im + stds*std_im

            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, extent=extent,
                            origin='lower', aspect='auto',
                            vmin = min, vmax = max)
        else:
            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, extent=extent,
                            origin='lower', aspect='auto')


        ax.set_title(kymo.get_plot_title())
        ax.set_ylabel('space (pixel)')
        ax.set_xlabel('time (s)')

        divider = make_axes_locatable(ax)
        fig.colorbar(im, cax = divider.append_axes('right', size='5%',
                    pad=0.05), orientation='vertical')

        plt.savefig(filename + kymo.get_file_title() + '.pdf', dpi=400,
                    bbox_inches='tight')
        plt.close()



# ===================================================== #
#  Butterworth filter  #
# ===================================================== #
def butter_bandpass(data, lowcut, highcut, fs, order):


    if lowcut==0 or lowcut==None:
        sos = signal.butter(order, highcut,  fs=fs, output='sos',
                            btype='lowpass')

    elif highcut==0 or highcut==None:
        sos = signal.butter(order, lowcut,  fs=fs, output='sos',
                            btype='highpass')

    else:
        sos = signal.butter(order, [lowcut, highcut], fs=fs, output='sos',
                            btype='bandpass')

    y = signal.sosfiltfilt(sos, data)
    return y



def butter_bandpass2d(data2d, lowcut, highcut, fs, order):
    filtered = data2d.copy()

    filtered = np.array([butter_bandpass(ts, lowcut, highcut, fs,
                        order) for ts in filtered])

    return filtered



def adaptive_butter_2d(data2d, fs, dom_freq, bandwidth, order,
                        min_lowcut=0.001):

    band_centers = np.unique(dom_freq)

    filtered = np.empty(( band_centers.size, *data2d.shape ))

    for i in range(band_centers.size):

        lowcut = np.max([band_centers[i]+bandwidth[0], min_lowcut ])
        highcut = band_centers[i]+bandwidth[1]

        filtered[i] = butter_bandpass2d(data2d, lowcut, highcut, fs, order)

    new = np.empty_like(data2d)
    time_window = data2d.shape[1]/dom_freq.size

    for i in range(dom_freq.size):

        choose_band, = np.where( band_centers==dom_freq[i] )
        start = int(i * time_window)
        stop = int((i+1) * time_window)

        new[:, start: stop] = filtered[choose_band[0], :, start: stop]

    return new

# ===================================================== #
#  time series  #
# ===================================================== #

def slid_aver(arr, window_size):
    pad_size = window_size
    expanded = np.pad(arr, pad_size , mode='reflect')
    new = np.convolve(expanded, np.ones(window_size)/window_size, mode='same')
    return new[pad_size:-pad_size]


def slid_aver2d(arr2d, window_size, axis=0):
    return np.apply_along_axis(slid_aver, axis, arr2d, *(window_size,))


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def plot_time_series(kymos, path=None, filename=None, window_size=10,
                    min_dist_s=70, spacial_avg=True):

    no_kymos = len(kymos)

    if np.any(path!=None):
        colspan=2
    else:
        colspan=3

    for i, kymo in zip(np.arange(no_kymos), kymos):
        ax = plt.subplot2grid((no_kymos, 3), (i,0), rowspan=1, colspan=colspan)

        if spacial_avg:
            time_series = slid_aver2d(kymo.kymo, window_size)[::window_size]
        else:
            time_series = kymo.kymo[::window_size]

        timestamps = np.arange(np.shape(time_series)[1]) *  kymos[0].t_scaling


        pixel_positions = np.arange(0, np.shape(time_series)[0] ) * window_size

        norm = mpl.colors.Normalize(vmin=pixel_positions.min(),
                                    vmax=pixel_positions.max())

        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        for t_s, p in zip(time_series, pixel_positions):
            min_p_dist = int(min_dist_s/ kymos[0].t_scaling)
            peak_inds, _ = signal.find_peaks(t_s, width=5, distance=min_p_dist)

            im = ax.plot(timestamps, t_s, c=cmap.to_rgba(p))

        ax.set_xlim(0, timestamps[-1])
        ax.set_ylabel(kymo.get_eq())
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.grid(True, axis='x', which='both')

        if not i==no_kymos-1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('time (s)')


    if np.any(path!=None):
        ax = plt.subplot2grid((no_kymos, 3), (0,2), rowspan=no_kymos, colspan=1)

        no_coords = path.shape[0]

        ax.scatter(path[:,1], path[:,0], c=range(no_coords), vmin=0,
                    vmax=no_coords, cmap='plasma')

        ax.axis('equal')
        ax.invert_yaxis()
        ax.axis('off')

    plt.subplots_adjust(hspace=0.5)
    plt.savefig(filename + 'time_series.pdf', dpi=400, bbox_inches='tight')
    plt.close()
    # plt.show()




# ===================================================== #
#  FFT  #
# ===================================================== #


def dominant_freq(kymo, lowcut=0.05, max_period=None):

    if max_period!=None:
        lowcut = 1./max_period

    f, Pxx_den_full = signal.periodogram(kymo.kymo, fs=1/kymo.t_scaling)
    Pxx_den = np.mean(Pxx_den_full, axis=0)

    Pxx_den = Pxx_den[f>lowcut]
    f = f[f>lowcut]

    dominant_freq = f[np.argmax(Pxx_den)]
    return dominant_freq



def power_spec(kymo, filename, min_frequency, lowcut=0.005, highcut=0.05,
                band=None, logscale=False):

    f, Pxx_den_full = signal.periodogram(kymo.kymo, fs=1./kymo.t_scaling,
                                scaling='spectrum')


    #### dominant frequency ->
    P = np.mean(Pxx_den_full, axis=0)
    P_cut = P[f>min_frequency]
    freq = f[f>min_frequency]
    dominant_freq = freq[np.argmax(P_cut)]
    #### <-


    # adjust for plotting..
    range = (f < highcut) * (f > lowcut)
    f          = f[range] * 1e3 ### switch to mHz
    Pxx_den_full    = Pxx_den_full[:,range]

    if logscale:
        Pxx_den_full = np.log(Pxx_den_full)


    Pxx_den         = np.mean(Pxx_den_full, axis=0)

    fig, axes = plt.subplots(2,1)
    ax = axes.ravel()

    # plotting...
    im = ax[0].imshow(Pxx_den_full, aspect='auto', cmap='plasma')

    cax = fig.add_axes([0.2, 0.95, 0.6, 0.05])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    if logscale:
        cbar.ax.set_title('log(power)', fontsize='medium')
    else:
        cbar.ax.set_title('power', fontsize='medium')

    # ax[0].set_title('log(power) ' + title)
    ax[0].set_xticks([])
    ax[0].set_ylabel('pixel')


    #
    ax[1].plot(f, Pxx_den)

    ax[1].axvline(x=dominant_freq*1e3, linewidth=1, color='red',
        label='dominant frequency: ' +
        str(np.around(dominant_freq*1e3, decimals=2)) + ' mHz')

    if np.all(band != None):
        mark1 = (dominant_freq + band[0]) * 1e3
        mark2 = (dominant_freq + band[1]) * 1e3
        ax[1].axvspan(mark1, mark2, facecolor='0.2',
                        alpha=0.3, label='frequency band')

    ax[1].legend()

    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].set_xlabel('frequency (mHz)')
    if logscale:
        ax[1].set_ylabel('log(power)')
    else:
        ax[1].set_ylabel('power')

    ax[1].set_xlim(f[0], f[-1])
    # fig.tight_layout()


    plt.savefig(filename + kymo.get_file_title() + '_power.pdf', dpi=200,
                bbox_inches='tight')
    plt.close()

    return dominant_freq



def spectro(kymo, filename, logscale=False, window_size_in_s=500, overlap=2,
             lowcut=0., highcut=0.04):

    window_size = int(window_size_in_s/kymo.t_scaling)

    window = signal.windows.blackmanharris(window_size)
    f, t, Sxx = signal.spectrogram(kymo.kymo, window = window,
                                    noverlap = window_size/overlap,
                                    fs=1./kymo.t_scaling, axis=-1)


    Sxx = np.mean(Sxx, axis=0)
    dom_freq = f[np.argmax(Sxx, axis=0)]


    #### plotting

    # crop
    range = (f < highcut) * (f > lowcut)
    f = f[range]
    Sxx    = Sxx[range]

    if logscale:
        Sxx=np.log(Sxx)
        label = 'log(power)'

    else:
        label = 'power'

    fig, ax = plt.subplots()

    ax.set_title('spectrogram ' + kymo.name)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('frequency (mHz)')


    extent =    [0, kymo.kymo.shape[1]*kymo.t_scaling,
                 0, f[-1]*1e3]

    im = ax.imshow(Sxx, cmap='plasma', extent=extent, origin='lower',
                    aspect='auto')

    divider = make_axes_locatable(ax)
    fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                 orientation='vertical', label=label)

    plt.savefig(filename + kymo.get_file_title() + '_spectro.pdf', dpi=200,
                bbox_inches='tight')
    plt.close()

    return dom_freq



# ===================================================== #
#  phase_shift2d   #
# ===================================================== #



# speeding up generic filter
def jit_filter_function(filter_function):
    jitted_function = numba.jit(filter_function, nopython=True)
    @cfunc(intc(CPointer(float64), intp, CPointer(float64), voidptr))
    def wrapped(values_ptr, len_values, result, data):
        values = carray(values_ptr, (len_values,), dtype=float64)
        result[0] = jitted_function(values)
        return 1
    return LowLevelCallable(wrapped.ctypes)



@jit_filter_function
def average_circular_dat(angles):
    """
    angles = set of angles in rad
    https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    """
    y = np.mean(np.sin(angles))
    x = np.mean(np.cos(angles))

    return np.arctan2(y, x)


def uni_filter_circular(arr, size_tuple):

    footprint=np.ones(size_tuple)
    return ndi.generic_filter(arr, average_circular_dat,
                                footprint=footprint, mode='reflect')




def phase_shift2d(phase1, phase2, filename, t_window_in_s=300):

    shift = phase2.kymo - phase1.kymo

    # map into  -pi, pi interval:
    shift = (shift + np.pi) % (2*np.pi) - np.pi

    #
    footprint = (20, int(t_window_in_s/phase1.t_scaling))
    avg_shift = uni_filter_circular(shift, footprint)


    ### plot ###
    # ........................................................................
    fig, axes = plt.subplots(2,1)
    ax = axes.ravel()
    plt.subplots_adjust(hspace=0.7)

    ax[0].set_title(r'phase difference $\Delta \phi\,$(rad)')
    ax[0].set_ylabel('space (pixel)')

    ax[1].set_title(r'phase difference $\langle \Delta \phi \rangle\,$(rad), footprint: '+ \
                    r'${0}$ pixels, '.format(footprint[0]) + \
                    r'${0}$ s, '.format(np.around(footprint[1]*phase1.t_scaling)))

    ax[1].set_xlabel('time (s)')
    ax[1].set_ylabel('space (pixel)')


    # limits
    extent =    [0, phase1.kymo.shape[1]*phase1.t_scaling,
                 0, phase1.kymo.shape[0]]

    # colorbar ticks
    ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
    labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']


    # twilight_r = cm.get_cmap('twilight_r', 100)
    # newcolors = twilight_r(np.linspace(0, 1, 100))
    # pink = np.array([248/256, 24/256, 148/256, 1])
    # newcolors[25:75, :] = pink
    # newcmp = ListedColormap(newcolors)

    im = ax[0].imshow(shift, extent=extent, origin='lower',
                    aspect='auto', cmap='twilight', vmin=-np.pi, vmax=np.pi)



    im = ax[1].imshow(avg_shift, extent=extent, origin='lower',
                    aspect='auto', cmap='twilight', vmin=-np.pi, vmax=np.pi)

    cbar = fig.colorbar(im, ax=ax.tolist(), orientation='vertical', ticks=ticks)
    cbar.ax.set_yticklabels(labels)

    plt.savefig(filename  + 'phase_shift_map.pdf', dpi=400, bbox_inches='tight')
    plt.close()

    return shift

# ===================================================== #
#  Phase dependence  #
# ===================================================== #

# binning #
def bin_acc_phase(phase, x, n_bins=10, norm=True):
    phase_f = phase.flatten()
    x_f     = x.flatten()

    x_mean, bin_edges,_   = stats.binned_statistic(phase_f, x_f,
                                                    statistic='mean',
                                                    bins=n_bins)
    x_std,_, _            = stats.binned_statistic(phase_f, x_f,
                                                    statistic=np.std,
                                                    bins=n_bins)

    if norm:
        x_max = np.max(x_mean)
        x_mean /= x_max
        x_std /= x_max
    return x_mean, bin_edges, x_std



def sin_func(phase, A, phi):
    return A * np.cos(phase - phi)


def phase_average(phase, kymos, filename, y_label, no_bins=15):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    phase_shift = []

    for k in kymos:
        bin_mean, bin_edges, bin_std = bin_acc_phase(phase, k.kymo,
                                                    n_bins=no_bins, norm=False)
        bin = bin_edges[:-1] + np.diff(bin_edges)/2


        phase_sample = np.linspace(bin[0], bin[-1], 100)
        popt, pcov   = curve_fit(sin_func, bin, bin_mean, p0=(-1,0))

        if k == kymos[0]:
            axis = ax1
            lab_text = k.get_eq(r'$(\phi)$')
        else:
            axis = ax2
            lab_text = k.get_eq(r'$(\phi)$') + ', shift: ' + \
                        str( np.around(popt[-1], decimals=2) )


        axis.plot(bin, bin_mean, label=lab_text, color=k.color)

        # plt.fill_between(bin, bin_mean - bin_std, bin_mean + bin_std,
        #                      alpha=0.2, color=k.color)

        # min = bin[np.argmin(bin_mean)]
        # plt.axvline(x=min, linewidth=1, color=k.color)

        axis.plot(phase_sample, sin_func(phase_sample, *popt), '--',
                    color=k.color)

        # y1 = sin_func(phase_sample, popt[0], popt[1] + pcov[1,1]**0.5)
        # y2 = sin_func(phase_sample, popt[0], popt[1] - pcov[1,1]**0.5)
        # plt.fill_between(phase_sample, y1, y2, color=k.color, alpha=0.15)

        axis.axvline(x=popt[-1], linewidth=1, color=k.color)
        phase_shift.append(popt[-1])


    ax1.set_ylabel('radius $(\mu$m$)$')
    ax2.set_ylabel(y_label)

    ax1.set_xlabel(r'phase $\phi$')
    # ax2.legend()
    ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
    labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center')

    plt.savefig(filename +'_'+ y_label + '_phase.pdf', dpi=400,
                bbox_inches='tight')
    plt.close()

    return phase_shift




# ===================================================== #
#  Correlation  #
# ===================================================== #

def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs


def correlation1d(kymo_ref, kymo_list, filename, title, upsample_t=1,
                        search_range_in_s=50):
    t_lims = int(search_range_in_s/(2*kymo_ref.t_scaling)) * upsample_t

    kymo_ref_zoom = ndi.zoom(kymo_ref.kymo, (1,upsample_t), order=5)
    delta_t = np.arange(-t_lims, t_lims+1)* kymo_ref.t_scaling / upsample_t

    fig, ax = plt.subplots()
    min_list = []
    for k in kymo_list:

        corr1d = match_template(ndi.zoom(k.kymo, (1,upsample_t), order=5),
                                kymo_ref_zoom[:,t_lims:-t_lims])[0]
        ind = np.unravel_index(np.argmin(corr1d, axis=None), corr1d.shape)
        min = ind[0] - int(corr1d.shape[0]/2)
        min_in_s = min * kymo_ref.t_scaling /upsample_t

        min_list.append(min_in_s)
        ax.plot(delta_t, corr1d, color=k.color)
        ax.axvline(x=min_in_s, linewidth=1, color=k.color,
                    label=r'$C$' + r'$[$' + k.get_eq('') +',' + \
                    kymo_ref.get_eq('')+ r'$](dt)$' + \
                    r', temporal shift $\Delta \mathcal{T}=$ ' \
                    + str(np.around(min_in_s, decimals=1)) + r'$\,$s')

    ax.set_xlabel(r'time lag $dt\,$(s)')
    ax.set_ylabel('norm. correlation')

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center')

    plt.savefig(filename + title + 'correlate1d.pdf', dpi=400,
                bbox_inches='tight')
    plt.close()
    return min_list


def correlation2d(kymo_a, kymo_b, filename, title,
                        upsample_t=1, upsample_x=1,
                        search_range_in_s=50, x_search=50):

    kymo1 = kymo_a.copy()
    kymo2 = kymo_b.copy()

    kymo1.kymo = ndi.zoom(kymo1.kymo, (upsample_x,upsample_t), order=5)
    kymo2.kymo = ndi.zoom(kymo2.kymo, (upsample_x,upsample_t), order=5)

    t_lims = int(search_range_in_s/(2*kymo1.t_scaling)) * upsample_t
    x_lims = np.min([upsample_x * x_search, int(kymo1.kymo.shape[0]/4.)])

    corr2d = match_template(kymo2.kymo,
                            kymo1.kymo[x_lims:-x_lims,t_lims:-t_lims])

    corr1d = corr2d[x_lims]

    ind = np.unravel_index(np.argmin(corr1d, axis=None), corr1d.shape)
    min = ind[0] - int(corr1d.shape[0]/2)
    min_in_s = min * kymo1.t_scaling / upsample_t

    #### IMSHOW PLOT ####
    fig, ax = plt.subplots()

    ext_t = t_lims * kymo1.t_scaling/ upsample_t
    im = ax.imshow(corr2d, aspect='auto',
                    extent=[-ext_t, ext_t, -x_lims, x_lims])

    ax.axvline(x=min_in_s, linewidth=1, color='grey',
                label='minimun along '+ r'$dx=0$')

    ax.set_xlabel(r'time lag $dt\,$(s)')
    ax.set_ylabel(r'space lag $dx\,$(pixels)')
    ax.set_title(r'$C[$' + kymo2.get_eq('') +',' + kymo1.get_eq('')+ r'$](dt, dx)$' )

    ax.legend()
    divider = make_axes_locatable(ax)
    fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                 orientation='vertical')

    fig.tight_layout()
    plt.savefig(filename + title + 'correlate2d.pdf', dpi=400,
                bbox_inches='tight')
    plt.close()

    return min_in_s


# ===================================================== #
#  Time/space dependend time shift #
# ===================================================== #
def time_shift2d(kymo_a, kymo_b, filename, upsample_t=1, window_size_in_s=200,
                    t_sampling_in_s = 10, x_sampling=20, search_range_in_s=100,
                    detrend=None):

    kymo1 = kymo_a.copy()
    kymo2 = kymo_b.copy()

    if detrend=='gauss':
        sec = 150
        sigma_gauss = sec / kymo1.t_scaling
        kymo1.kymo = gauss_detrend(kymo1.kymo, sigma_gauss)
        kymo2.kymo = gauss_detrend(kymo2.kymo, sigma_gauss)

    window_size = np.around(window_size_in_s/kymo1.t_scaling).astype(int)
    search_extend = np.around(search_range_in_s/(2.*kymo1.t_scaling)).astype(int)
    t_sampling = np.around(t_sampling_in_s/ kymo1.t_scaling ).astype(int)
    x_sampling = np.min([x_sampling, kymo1.kymo.shape[0]])

    s_min = search_extend * upsample_t
    s_max = -search_extend * upsample_t

    window  = (x_sampling, (window_size + 2*search_extend)*upsample_t )
    step    = (x_sampling, int(t_sampling * upsample_t) )

    kymo1.kymo = ndi.zoom(kymo1.kymo, (1,upsample_t), order=5)
    kymo2.kymo = ndi.zoom(kymo2.kymo, (1,upsample_t), order=5)

    a_windows = view_as_windows(kymo1.kymo, window_shape=window, step=step)
    b_windows = view_as_windows(kymo2.kymo, window_shape=window, step=step)
    shape = np.shape(b_windows)

    mins = np.zeros(shape[:2])

    for i in range(shape[0]):
        for j in range(shape[1]):
            # crop window to orig. window_size
            a = a_windows[i,j][:, s_min:s_max]
            b = b_windows[i,j]

            corr1d = match_template(b,a)[0]

            ind = np.unravel_index(np.argmin(corr1d, axis=None), corr1d.shape)
            mins[i,j] = ind[0] - int(corr1d.shape[0]/2)
            # if np.abs(mins[i,j]) >= s_min-1:
            # show_im(a)
            # show_im(b)
            #
            # plt.plot(corr1d)
            # plt.show()
            # plt.close()


    time_shift_map = mins * kymo1.t_scaling/ upsample_t

    # plotting...
    fig, ax = plt.subplots()

    ax.set_title(r'time shift $\Delta \mathcal{T}\,$(s)')
    ax.set_ylabel('space (pixel)')
    ax.set_xlabel('time (s)')

    # limits
    extent =    [0, kymo_a.kymo.shape[1]*kymo_a.t_scaling,
                 0, kymo_a.kymo.shape[0]]

    vminmax = np.max(np.abs([np.min(time_shift_map), np.max(time_shift_map)]))

    im = ax.imshow(time_shift_map, extent=extent, origin='lower',
                    aspect='auto', cmap='twilight_shifted',
                    vmin=-search_range_in_s/2.,
                    vmax=+search_range_in_s/2.)


    # colorbar ticks
    ticks = np.arange(-100, 100, 10)
    labels = ticks.astype(str)

    mean = np.mean(time_shift_map)
    hide_id = np.argmin(np.abs(ticks - mean))
    labels[hide_id] = ''

    ticks = np.append(ticks, mean)
    labels = np.append(labels,
                r'mean: ${0}$'.format(str(np.around(mean, decimals=1))))


    # colorbar
    divider = make_axes_locatable(ax)
    cbar = fig.colorbar(im, cax = divider.append_axes('right', size='5%',
                        pad=0.05), orientation='vertical', ticks=ticks)
    cbar.ax.set_yticklabels(labels)


    plt.savefig(filename + kymo1.state + 'time_shift_map.pdf', dpi=400,
                bbox_inches='tight')
    plt.close()

    return time_shift_map





# ===================================================== #
#  analyt. flow calc. #
# ===================================================== #

def estimate_flow(radii, frame_int):
    radii_b =  bandpass(radii, frame_int, highcut=0.015)
    integ = np.gradient(radii_b, axis=1) * radii_b

    flow = np.cumsum(integ, axis=0)
    flow /= np.max(flow)

    plt.plot(flow[1])
    plt.plot(flow[-1])

    # plt.show()

    show_im(flow[:, 500:-100])

    return
