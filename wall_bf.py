"""

"""

from analysis.network_analysis import *
from analysis.data_sets import *
import skimage.morphology as morph
import os
from analysis.plotting import *

import matplotlib.animation as animation
import matplotlib.patches as mpatches
import matplotlib.cm as mc
import matplotlib.colors as mcol
from skimage import feature

import multiprocessing
import itertools



def smooth_area(im, sigma, thresh):
    area = ndi.gaussian_filter(im, sigma=sigma)
    area = np.where(area > thresh, True, False)
    area = morph.remove_small_holes(area, area_threshold=200, connectivity=2)
    area = morph.remove_small_objects(area, min_size=1000, connectivity=2)
    return area.astype(float)


def detect_wall(i, data_sets, ref_range, plot_it=False):

    set = data_sets[i]
    print("Detect wall: ", set.file_dat)

    bf      = invert_bf(read_file(set.file_raw))
    bf_ims  = [invert_bf(read_file(ref.file_raw)) for ref in data_sets[i-ref_range:i+ref_range]]

    mask = create_mask(bf, set.sigma, set.threshold, set.halo_sig)
    ectoplasm = smooth_area(mask, 40, 0.2)

    skeleton = extract_skeleton(ectoplasm, method='medial_axis', branch_thresh=150)
    r_ecto = np.sum(extract_radii(ectoplasm, skeleton))/np.sum(skeleton)


    diff = np.var(bf_ims, axis=0)
    motion = np.where(diff > np.mean(diff), 1., 0.)

    endoplasm = smooth_area(motion, 40, 0.3)
    r_endo = np.sum(extract_radii(endoplasm, skeleton))/np.sum(skeleton)

    wall_map = 2 * ectoplasm + endoplasm

    np.savez_compressed(set.file_dat,   wall_map = wall_map,
                                        r_ecto = r_ecto,
                                        r_endo = r_endo)

    if plot_it:
        plt.imshow(bf/np.max(bf), cmap='gray')
        plt.imshow(endoplasm + ectoplasm*2, alpha=0.1)
        plt.imshow(thick_skeleton(skeleton))
        plt.savefig(set.file_plot + 'wall.png', dpi=100)
        plt.close()



def collect_data_from_npz(data_sets, indices, key):
    values = []
    for i in indices:
        values.append(np.load(data_sets[i].file_dat + '.npz')[key])

    return np.array(values)



def anim_npz(data_sets, frame_step):

    data_sets = data_sets[::frame_step]
    set = data_sets[0]

    timestamps = np.arange(len(data_sets)) * frame_step * data_sets[0].frame_int



    fig, ax=plt.subplots()

    # def cmap, get colors, def labels
    ref_im = np.load(data_sets[0].file_dat + '.npz')['wall_map']
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
    ani.save(data_sets[0].file_plot_set + 'wall_bf.mp4')


#########################################
############### wall_bf #################
#########################################

def main(): ## python3 wall_bf.py <keyword>

    set_keyword = os.sys.argv[1].strip()

    data_set_range = range(data(set_keyword).first, data(set_keyword).last)
    # data_set_range = range(1, 13)


    data_sets = [data(set_keyword, i, method='wall', color='bf') for i in data_set_range]

    frame_step = 10
    ref_range = 5


    ########### ANALYSE #############
    index_ar = range(ref_range, len(data_sets)-ref_range)

    num_threads = multiprocessing.cpu_count()
    print(num_threads)
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



    ########### PLOTTING #############
    print("Plotting...")
    timestamps = np.arange(len(r_ectos)) * data_sets[0].frame_int

    norm = np.max(r_ectos)
    plt.plot(timestamps, r_ectos/norm, label = 'tube')
    plt.plot(timestamps, r_endos/norm, label = 'endoplasm')
    plt.plot(timestamps, r_endos/r_ectos, label = 'ratio endo/ecto')

    plt.xlabel('time (s)')
    plt.ylabel('radius (a.u.)')
    plt.legend()

    plt.savefig(data_sets[0].file_plot_set + 'wall_bf.pdf')
    plt.close()


    ########### ANIMATE #############
    anim_npz(data_sets[ref_range:-ref_range], frame_step)



if __name__ == "__main__":
    main()
