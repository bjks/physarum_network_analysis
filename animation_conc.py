import matplotlib.animation as animation
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
import matplotlib.cm as mc

import os


def anim(ims, set, frame_step, keyword):

    if 'raw' in keyword:
        std_im, mean_im = np.nanstd(ims[0]), np.nanmean(ims[0])
        min, max = mean_im - std_im, mean_im + std_im
    else:
        min = np.nanmin(ims[0])
        max = np.nanmax(ims[0])
        print(min, max)

    timestamps = np.arange(len(ims)) * frame_step * set.frame_int

    fig, ax=plt.subplots()
    cmap = mc.get_cmap('Spectral_r')

    images = []
    for im, t in zip(ims, timestamps):
        title = ax.text(0.5,1.03,"time: {0:.2f}s".format(t),
                size=plt.rcParams["axes.titlesize"],
                ha="center", transform=ax.transAxes)

        temp_im = ax.imshow(im, animated=True, cmap=cmap,
                            vmin = min, vmax = max)
        images.append([temp_im, title])

    plt.colorbar(images[0][0])

    ani = animation.ArtistAnimation(fig, images, interval=200, blit=True)
    if not os.path.exists(set.path_plots):
        os.mkdir(set.path_plots)
    ani.save(set.file_plot_set + keyword + '.mp4', dpi=200)

#################################################################

def main():

    set_keyword = os.sys.argv[1]
    keyword     = os.sys.argv[2]
    color       = os.sys.argv[3]
    step        = int(os.sys.argv[4])

    method = 'inter_mean'

    first, last = data(set_keyword).first, data(set_keyword).last
    t_arr = np.arange(first,last, step)

    data_sets = [data(set_keyword, i, method, color=color) for i in t_arr]

    #################################################################

    map = []
    print('>> Reading...')

    for set, t in zip(data_sets, t_arr):
        print('>>', t)

        if keyword == 'raw':
            green_clean = np.load(set.file_dat + '.npz')['green_clean']
            texas_clean = np.load(set.file_dat + '.npz')['texas_clean']
            image = calc_ratio(green_clean, texas_clean)


        elif keyword == 'raw_net':
            image = np.load(set.file_dat + '.npz')['network_clean']


        elif keyword == 'raw_norm':
            rd          = np.load(set.file_dat + '.npz')['rel_dist']
            mask        = np.load(set.file_dat + '.npz')['mask']
            local_radii = np.load(set.file_dat + '.npz')['local_radii']
            skeleton    = np.load(set.file_dat + '.npz')['skeleton']

            r  = tube_radius_at_point(mask, skeleton, local_radii)
            height = np.where(rd<1, r * np.sqrt(1 - rd**2), 0)
            if color == 'green':
                network = np.load(set.file_dat + '.npz')['network_clean']
            else:
                network = np.load(set.file_dat + '.npz')['green_clean']

            image = np.where(height>0, calc_ratio(network,height), 0)


        elif keyword == 'conc_rat':
            concentration   = np.load(set.file_dat + '.npz')['concentration']
            mask            = np.load(set.file_dat + '.npz')['mask']
            skeleton        = np.load(set.file_dat + '.npz')['skeleton']

            image = tube_radius_at_point(mask, skeleton, concentration)


        elif keyword == 'norm_conc':
            concentration   = np.load(set.file_dat + '.npz')['concentration']
            local_radii     = np.load(set.file_dat + '.npz')['local_radii']
            mask            = np.load(set.file_dat + '.npz')['mask']

            norm_conc = calc_ratio(concentration, local_radii)

            image = tube_radius_at_point(mask, skeleton, norm_conc)


        elif keyword == 'skeleton':
            mask            = np.load(set.file_dat + '.npz')['mask'].astype(int)
            skeleton        = np.load(set.file_dat + '.npz')['skeleton'].astype(int)

            selem = morph.disk(5)
            thick_skeleton = morph.dilation(skeleton, selem)

            image = thick_skeleton + mask



        image = np.where(image == 0, np.nan, image)
        map.append(image)

    print('>> Create animation...')
    anim(map, data_sets[0], step, keyword)


if __name__ == "__main__":
    main()
