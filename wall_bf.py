"""

"""

from analysis.network_analysis import *
from analysis.data_sets import *
import skimage.morphology as morph
import os
import sys
from analysis.plotting import *

import matplotlib.animation as animation
import matplotlib.patches as mpatches
import matplotlib.cm as mc
import matplotlib.colors as mcol
from skimage import feature
from skimage import transform

from skimage.util.shape import view_as_windows

import multiprocessing
import itertools

import time



def collect_data_from_npz(data_sets, indices, key):
    values = []
    for i in indices:
        values.append(np.load(data_sets[i].file_dat + '.npz')[key])

    return np.array(values)



def sliding_avg(arr, window_size, mode='valid'):
    return np.convolve(arr, np.ones(window_size)/window_size, mode=mode)



def remove_outlier(arr, stds):
    range = stds * np.nanstd(arr)
    mean  = np.nanmean(arr)

    return np.where((arr < mean+range) * (arr > mean-range), arr, np.nan)



def bring_to_shape_of(ar, ar_with_desired_shape):
    return transform.resize(ar, np.shape(ar_with_desired_shape))



def signed_avg(a, b, for_sign=1):
    return np.sqrt( np.dot(a,a) + np.dot(b,b) ) * np.sign(for_sign)



def smooth_area(im, sigma):
    area = ndi.gaussian_filter(im, sigma=sigma)

    thresh = threshold_otsu(area)

    area = np.where(area > thresh, True, False)

    ##### Some corrections ...: #######
    area = morph.remove_small_holes(area, area_threshold=int(1e4), connectivity=2)
    # area = morph.remove_small_objects(area, min_size=2000, connectivity=2)

    area = extract_network(area, 1)

    return area.astype(float)



def flow_quantification(frame_a, frame_b, mask, upsample=2, window_size=40, sampling=20):

    window  = (window_size * upsample,) * 2
    step    = sampling * upsample

    if upsample > 1:
        mask    = ndi.zoom(mask.astype(bool), upsample, order=0)         # order 0!!
        frame_a = ndi.zoom(frame_a, upsample, order=5)
        frame_b = ndi.zoom(frame_b, upsample, order=5)

    frame_a = np.multiply(frame_a, mask)
    frame_b = np.multiply(frame_b, mask)

    shape = np.shape(view_as_windows(frame_a, window_shape= window, step=step))

    flow_field_x = np.zeros(shape[:2])
    flow_field_y = np.zeros(shape[:2])

    for i in range(shape[0]):
        for j in range(shape[1]):
            a = view_as_windows(frame_a, window_shape= window, step=step)[i,j]
            b = view_as_windows(frame_b, window_shape= window, step=step)[i,j]

            image_product = np.fft.fft2(a) * np.fft.fft2(b).conj()
            cc_image = np.fft.ifft2(image_product)

            corr = np.fft.fftshift(cc_image).real

            if np.max(corr)>0:
                shift = np.unravel_index(np.argmax(corr, axis=None), corr.shape)

                flow_field_x[i,j] = shift[1]  - window[1]/2
                flow_field_y[i,j] = shift[0]  - window[0]/2


    flow_field_x/=upsample
    flow_field_y/=upsample

    flow_x = np.nanmean(flow_field_x)
    flow_y = np.nanmean(flow_field_y)


    avg_flow = signed_avg(flow_x, flow_y, for_sign = flow_x)

    return avg_flow, flow_field_x, flow_field_y



def detect_wall(i, data_sets, ref_range, plot_it=False):

    set = data_sets[i]
    print("Detect wall: ", set.file_dat)

    ### reference image and image for flow quantification
    bf      = invert_bf(read_file(set.file_raw))
    bf_next = invert_bf(read_file(data_sets[i+1].file_raw))


    ### calc ectoplasm based on simple mask
    mask = create_mask(bf, set.sigma, set.threshold, set.halo_sig)
    ectoplasm = smooth_area(mask, 30)

    ### calc radii based on medial axis in ectoplasm map
    skel_ecto = extract_skeleton(ectoplasm, method='skeletonize', branch_thresh=250)
    r_ecto = np.sum(extract_radii(ectoplasm, skel_ecto))/np.sum(skel_ecto)


    ### calc endoplasm based on temporal variance in each point
    diff = np.var([invert_bf(read_file(ref.file_raw)) for ref in data_sets[i-ref_range:i+ref_range]], axis=0)
    motion = np.where(diff > np.mean(diff), 1., 0.)

    endoplasm = smooth_area(motion, 30)

    ### calc radii ...
    skel_endo = extract_skeleton(endoplasm, method='skeletonize', branch_thresh=250)
    r_endo = np.sum(extract_radii(endoplasm, skel_endo))/np.sum(skel_endo)

    wall_map = 2 * ectoplasm + endoplasm

    avg_flow, flow_field_x, flow_field_y = flow_quantification(bf, bf_next, ectoplasm)


    np.savez_compressed(set.file_dat,   wall_map        = wall_map,
                                        skel_ecto       = skel_ecto,
                                        skel_endo       = skel_endo,
                                        r_ecto          = r_ecto,
                                        r_endo          = r_endo,
                                        avg_flow        = avg_flow,
                                        flow_field_x    = flow_field_x,
                                        flow_field_y    = flow_field_y)




################################################################################
def analyse_sets(index_ar, data_sets, ref_range, debug_mode=False):
    if debug_mode:
        for i in index_ar:
            detect_wall(i, data_sets, ref_range)

    else:
        num_threads = multiprocessing.cpu_count()
        print("Number of detected cores: ", num_threads)

        p = multiprocessing.Pool(num_threads)
        p.starmap(detect_wall, zip( index_ar,
                                    itertools.repeat(data_sets),
                                    itertools.repeat(ref_range)))
        p.close()
        p.join()


    ########### COLLECT #############
    print("Collect data from npz files ...")
    r_endos     = collect_data_from_npz(data_sets, index_ar, 'r_endo')
    r_ectos     = collect_data_from_npz(data_sets, index_ar, 'r_ecto')
    avg_flow    = collect_data_from_npz(data_sets, index_ar, 'avg_flow')


    ########## SAVE SUMMARY #########
    np.savez_compressed(data_sets[0].file_dat_set,  index_ar    = index_ar,
                                                    ref_range   = ref_range,
                                                    r_ectos     = r_ectos,
                                                    r_endos     = r_endos,
                                                    avg_flow    = avg_flow)




################################################################################
def plot_dynamics(data_sets, window_sli_avg):
    r_ectos = np.load(data_sets[0].file_dat_set + '.npz')['r_ectos']
    r_endos = np.load(data_sets[0].file_dat_set + '.npz')['r_endos']
    avg_flow = np.load(data_sets[0].file_dat_set + '.npz')['avg_flow']

    r_ectos = sliding_avg(r_ectos, window_sli_avg)

    r_endos = remove_outlier(r_endos, 1)
    r_endos = sliding_avg(r_endos, window_sli_avg)
    avg_flow = sliding_avg(avg_flow, window_sli_avg)



    print("Plotting...")
    timestamps = np.arange(len(r_ectos)) * data_sets[0].frame_int

    ########## WALL #########
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('normalized radii (a.u.)')

    ax2.set_ylim([0.44, 0.75])
    ax2.set_ylabel('ratio endo/ecto')


    lns1 = ax1.plot(timestamps, r_ectos/np.nanmax(r_ectos),
                    label = 'tube radius', color='mediumblue')
    lns2 = ax1.plot(timestamps, r_endos/np.nanmax(r_endos),
                    label = 'endoplasm radius', color='lightskyblue')
    lns3 = ax2.plot(timestamps, r_endos/r_ectos,
                    label = 'ratio endo/ecto', color='orange')
    # ax1.grid()

    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs)

    plt.savefig(data_sets[0].file_plot_set + 'wall.pdf')
    plt.close()



    ########## FLOW #########
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('normalized flow (a.u.)')
    ax2.set_ylabel('ratio endo/ecto')

    lns1 = ax1.plot(timestamps, avg_flow/np.nanmax(avg_flow),
                    label = 'flow', color='green')
    lns2 = ax2.plot(timestamps, r_endos/r_ectos,
                    label = 'ratio endo/ecto', color='orange')


    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs)

    plt.savefig(data_sets[0].file_plot_set + 'flow.pdf')
    plt.close()



################################################################################
def plot_samples(data_sets, index):

    for i in index:
        set = data_sets[i]
        dat = np.load(set.file_dat + '.npz')

        bf      = invert_bf(read_file(set.file_raw))

        plt.imshow(bf/np.max(bf), cmap='gray')
        plt.imshow(dat['wall_map'], alpha=0.1)
        plt.imshow(thick_skeleton(dat['skel_ecto']))

        plt.savefig(set.file_plot + 'wall.png', dpi=200)
        plt.close()

        abs_flow = bring_to_shape_of(np.sqrt(   dat['flow_field_x']**2 + \
                                                dat['flow_field_y']**2      ), \
                                                dat['wall_map']             )

        plt.imshow(bf/np.max(bf), cmap='gray')
        plt.imshow(abs_flow, alpha=0.2)
        plt.colorbar()
        plt.savefig(set.file_plot + 'flow.png', dpi=200)
        plt.close()




################################################################################
def anim_npz(data_sets, frame_step):

    data_sets = data_sets[::frame_step]
    set = data_sets[0]

    timestamps = np.arange(len(data_sets)) * frame_step * data_sets[0].frame_int

    fig, ax=plt.subplots()

    # def cmap, get colors, def labels
    ref_im = np.load(data_sets[0].file_dat + '.npz')['wall_map']
    # values = [0,1,3]
    values = np.unique(ref_im.ravel())

    cmap = mc.get_cmap('Blues')
    norm = mcol.Normalize(vmin=np.min(values), vmax=np.max(values))
    colors = [cmap(norm(value)) for value in values]
    labels = ['', 'ectoplasm', 'endoplasm']
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(1,len(values))]
    ax.legend(handles=patches)


    images = []
    for d, t in zip(data_sets, timestamps):
        title = ax.text(0.5,1.03,"time: {0:.2f}s".format(t),
                size=plt.rcParams["axes.titlesize"],
                ha="center", transform=ax.transAxes)

        temp_map = np.load(d.file_dat + '.npz')['wall_map']

        temp_im = ax.imshow(temp_map, animated=True, cmap=cmap)
        images.append([temp_im, title])

    ani = animation.ArtistAnimation(fig, images, interval=200, blit=True)
    if not os.path.exists(data_sets[0].path_plots):
        os.mkdir(data_sets[0].path_plots)
    ani.save(data_sets[0].file_plot_set + '.mp4')



################################################################################
############################## MAIN ############################################
################################################################################
def main(): ## python3 wall_bf.py <keyword>


#########################################
#########################################
    set_keyword     = os.sys.argv[1].strip()

    ref_range       = 20
    window_sli_avg  = 20


    if sys.platform.startswith('linux'):
        analyse_bool        = True
        plot_dynamics_bool  = True
        animate_bool        = True
        plot_step           = 100

        frame_step = 15
        data_set_range = range(data(set_keyword).first, data(set_keyword).last)


    else:
        analyse_bool        = False
        plot_dynamics_bool  = True
        animate_bool        = False
        plot_step           = 1

        frame_step = 1
        data_set_range = range(1, 44)

#########################################
#########################################


    data_sets = [data(set_keyword, i, method='wall', color='bf') for i in data_set_range]
    index_ar = range(ref_range, len(data_sets)-ref_range)

    if analyse_bool:
        analyse_sets(index_ar, data_sets, ref_range)


    if plot_dynamics_bool:
        plot_dynamics(data_sets, window_sli_avg)


    if plot_step != None:
        inds_plot = np.arange(ref_range, len(data_sets)-ref_range, plot_step)
        plot_samples(data_sets, inds_plot)


    if animate_bool:
        anim_npz(data_sets[ref_range:-ref_range], frame_step)



if __name__ == "__main__":
    main()
