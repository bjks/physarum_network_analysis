import matplotlib.animation as animation
import matplotlib.cm as mc
import os

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *



def anim_arr(ims, set, frame_step, keyword):

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
#################################################################


def anim_npz(data_sets, frame_step, keyword):
    frame_step = np.max([frame_step,1])
    data_sets = data_sets[::frame_step]
    set = data_sets[0]

    timestamps = np.arange(len(data_sets)) * frame_step * data_sets[0].frame_int
    fig, ax=plt.subplots()
    colorbar = True


    images = []
    for set, t in zip(data_sets, timestamps):
        print(set.file_dat)

        title = ax.text(0.5,1.03,"time: {0:.2f}s".format(t),
                size=plt.rcParams["axes.titlesize"],
                ha="center", transform=ax.transAxes)


        if keyword == 'raw':
            green_clean = np.load(set.file_dat + '.npz')['green_clean']
            texas_clean = np.load(set.file_dat + '.npz')['texas_clean']

            image = calc_ratio(green_clean, texas_clean)
            temp_map = np.where(image == 0, np.nan, image)

            if set == data_sets[0]:
                std_im, mean_im = np.nanstd(temp_map), np.nanmean(temp_map)

                print(mean_im, std_im)
                min, max = mean_im - 2*std_im, mean_im + 2*std_im

            cmap = mc.get_cmap('Spectral_r')

            temp_im = ax.imshow(temp_map, animated=True, cmap=cmap,
                                vmin = min, vmax = max)


        elif keyword == 'raw_net':
            image = np.load(set.file_dat + '.npz')['network_clean']

            if set == data_sets[0]:
                std_im, mean_im = np.nanstd(image), np.nanmean(image)
                min, max = mean_im - std_im, mean_im + std_im

            cmap = mc.get_cmap('Spectral_r')

            temp_map = np.where(image == 0, np.nan, image)
            temp_im = ax.imshow(temp_map, animated=True, cmap=cmap,
                                vmin = min, vmax = max)

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

            if set == data_sets[0]:
                min, max = np.nanmin(image), np.nanmax(image)

            cmap = mc.get_cmap('Spectral_r')

            temp_map = np.where(image == 0, np.nan, image)
            temp_im = ax.imshow(temp_map, animated=True, cmap=cmap,
                                vmin = min, vmax = max)


        elif keyword == 'conc_rat':
            concentration   = np.load(set.file_dat + '.npz')['concentration']
            mask            = np.load(set.file_dat + '.npz')['mask']
            skeleton        = np.load(set.file_dat + '.npz')['skeleton']

            image = tube_radius_at_point(mask, skeleton, concentration)
            if set == data_sets[0]:
                std_im, mean_im = np.nanstd(image), np.nanmean(image)
                min, max = mean_im - std_im, mean_im + std_im

            if set == data_sets[0]:
                min, max = np.nanmin(image), np.nanmax(image)

            cmap = mc.get_cmap('Spectral_r')

            temp_map = np.where(image == 0, np.nan, image)
            temp_im = ax.imshow(temp_map, animated=True, cmap=cmap,
                                vmin = min, vmax = max)



        elif keyword == 'norm_conc':
            concentration   = np.load(set.file_dat + '.npz')['concentration']
            local_radii     = np.load(set.file_dat + '.npz')['local_radii']
            mask            = np.load(set.file_dat + '.npz')['mask']

            norm_conc = calc_ratio(concentration, local_radii)

            image = tube_radius_at_point(mask, skeleton, norm_conc)
            if set == data_sets[0]:
                min, max = np.nanmin(image), np.nanmax(image)

            cmap = mc.get_cmap('Spectral_r')

            temp_map = np.where(image == 0, np.nan, image)
            temp_im = ax.imshow(temp_map, animated=True, cmap=cmap,
                                vmin = min, vmax = max)


        elif keyword == 'skeleton':
            mask        = np.load(set.file_dat + '.npz')['mask'].astype(int)
            skeleton    = np.load(set.file_dat + '.npz')['skeleton'].astype(int)

            selem = morph.disk(5)
            thick_skeleton = morph.dilation(skeleton, selem)

            image = mask + thick_skeleton
            cmap = mc.get_cmap('Blues')
            temp_im = ax.imshow(image, animated=True, cmap=cmap)
            colorbar = False


        elif keyword == 'mask':
            mask = np.load(set.file_dat + '.npz')['mask'].astype(int)
            cmap = mc.get_cmap('Blues')
            temp_im = ax.imshow(mask, animated=True, cmap=cmap)
            colorbar = False


        elif keyword == 'flow':
            x = np.load(set.file_dat + '.npz')['flow_x']
            y = np.load(set.file_dat + '.npz')['flow_y']

            temp_im = ax.quiver(x,y, scale=1, scale_units='x')
            colorbar = False

        images.append([temp_im, title])

    if colorbar:
        fig.colorbar(temp_im)
    ani = animation.ArtistAnimation(fig, images, interval=200, blit=True)
    if not os.path.exists(data_sets[0].path_plots):
        os.mkdir(data_sets[0].path_plots)
    ani.save(data_sets[0].file_plot_set + '_' + keyword + '_.mp4', dpi=200)



#################################################################
#################################################################
#################################################################
#################################################################

def main():

    set_keyword = os.sys.argv[1]
    keyword     = os.sys.argv[2]
    color       = os.sys.argv[3]
    no          = int(os.sys.argv[4])

    if no > 300:
        step = 1
        data_set_range = range(data(set_keyword).first, no)

    else:
        step = no
        data_set_range = range(data(set_keyword).first, data(set_keyword).last)


    if keyword=='flow':
        method = 'piv'
    else:
        method = 'inter_mean'

    # data_set_range = range(1, 13)

    data_sets = [data(set_keyword, i, method=method, color=color) for i in data_set_range]

    ####### ANIMATE #######
    anim_npz(data_sets, step, keyword)

if __name__ == "__main__":
    main()
