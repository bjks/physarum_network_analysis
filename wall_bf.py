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


def smooth_area(im, sigma, thresh):
    area = ndi.gaussian_filter(im, sigma=sigma)
    area = np.where(area > thresh, True, False)
    area = morph.remove_small_holes(area, area_threshold=200, connectivity=2)
    area = morph.remove_small_objects(area, min_size=1000, connectivity=2)
    return area.astype(float)


def detect_wall(ref_sets, set, plot_it=False):

    bf      = invert_bf(read_file(set.file_raw))
    bf_ims  = [invert_bf(read_file(ref.file_raw)) for ref in ref_sets]

    mask = create_mask(bf, set.sigma, set.threshold, set.halo_sig)
    ectoplasm = smooth_area(mask, 40, 0.2)

    skeleton = extract_skeleton(ectoplasm, method='medial_axis', branch_thresh=150)
    r_ecto = np.sum(extract_radii(ectoplasm, skeleton))/np.sum(skeleton)


    diff = np.var(bf_ims, axis=0)
    motion = np.where(diff > np.mean(diff), 1., 0.)

    endoplasm = smooth_area(motion, 40, 0.3)
    r_endo = np.sum(extract_radii(endoplasm, skeleton))/np.sum(skeleton)
    if plot_it:
        plt.imshow(bf/np.max(bf), cmap='gray')
        plt.imshow(endoplasm + ectoplasm*2, alpha=0.1)
        plt.imshow(thick_skeleton(skeleton))
        plt.savefig(set.file_plot + 'wall.png', dpi=100)
        plt.close()
    return 2 * ectoplasm + endoplasm, r_ecto, r_endo



def anim(ims, set, frame_step):

    timestamps = np.arange(len(ims)) * frame_step * set.frame_int

    fig, ax=plt.subplots()

    # def cmap, get colors, def labels
    values = np.unique(ims[0].ravel())
    cmap = mc.get_cmap('Blues')
    norm = mcol.Normalize(vmin=np.min(values), vmax=np.max(values))
    colors = [cmap(norm(value)) for value in values]
    labels = ['', 'ectoplasm', 'endoplasm']
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(1,len(values))]
    ax.legend(handles=patches)


    images = []
    for im, t in zip(ims, timestamps):
        title = ax.text(0.5,1.03,"time: {0:.2f}s".format(t),
                size=plt.rcParams["axes.titlesize"],
                ha="center", transform=ax.transAxes)

        temp_im = ax.imshow(im, animated=True, cmap=cmap)
        images.append([temp_im, title])

    ani = animation.ArtistAnimation(fig, images, interval=200, blit=True)
    if not os.path.exists(set.path_plots):
        os.mkdir(set.path_plots)
    ani.save(set.file_plot_set + 'wall_bf.mp4')



#########################################
############### wall_bf #################
#########################################

def main(): ## python3 wall_bf.py <keyword>

    set_keyword = os.sys.argv[1].strip()
    data_sets = [data(set_keyword, i, color='bf') for i in range(data(set_keyword).first, data(set_keyword).last)]
    # data_sets = [data(set_keyword, i, color='bf') for i in range(1, 13)]

    frame_step = 20
    ref_range = 5


    wall_maps = []
    r_endos = []
    r_ectos = []

    for i in range(ref_range, len(data_sets)-ref_range):
        print("Wall", i, set_keyword)
        map, r_ecto, r_endo = detect_wall(data_sets[i-ref_range:i+ref_range], data_sets[i])
        if i%frame_step==0:
            wall_maps.append(map)
        r_endos.append(r_endo)
        r_ectos.append(r_ecto)



    timestamps = np.arange(len(r_ectos)) * data_sets[0].frame_int

    norm = np.max(r_ectos)
    plt.plot(timestamps, r_ectos/norm, label = 'tube')
    plt.plot(timestamps, r_endos/norm, label = 'endoplasm')

    plt.xlabel('time (s)')
    plt.ylabel('radius (a.u.)')
    plt.legend()

    plt.savefig(data_sets[0].file_plot_set + 'wall_bf.pdf')
    plt.close()

    anim(wall_maps, data_sets[0], frame_step)



if __name__ == "__main__":
    main()
