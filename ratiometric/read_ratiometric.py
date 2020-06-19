import matplotlib.pyplot as plt
import matplotlib as mpl

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *
from ratiometric import *
from skeleton.skeleton import *
from analysis.tools import *


from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar # scalebar
plt.rcParams["font.family"] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Avenir', 'CMU Sans Serif']



def processing_visualization(file_dat, file_plot, analyse_flow,show=False, pretty=False):
    print("Plotting ", file_dat)
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(2,3, sharex=True, sharey=True)
    ax = axes.ravel()


    green           = data['green_clean']
    green = np.ma.masked_where(green == 0, green)

    texas           = data['texas_clean']
    texas = np.ma.masked_where(texas == 0, texas)

    radius = np.ma.masked_where(data['radii_map'] == 0, data['radii_map'])
    ax[0].imshow(radius)
    ax[0].set_title('mapped radius')

    ax[1].imshow(green, cmap='Greens') #, vmin = min, vmax = max)
    ax[1].set_title('green channel')

    ax[2].imshow(texas, cmap='Reds')
    ax[2].set_title('red channel')

    if analyse_flow:
        abs_flow = np.sqrt(data['flow_x']**2 +  data['flow_y']**2)
        ax[3].imshow(thick_skeleton(abs_flow))
        ax[3].set_title('flow speed')

        ax[0].imshow(thick_skeleton(data['skeleton']), cmap='Reds_r')

    else:
        ax[3].set_title('skeleton')
        ax[3].imshow(thick_skeleton(data['skeleton']), cmap='Blues_r')



    ax[4].imshow(thick_skeleton(data['c_green']), cmap='Greens')
    ax[4].set_title('projected green')

    ax[5].imshow(thick_skeleton(data['c_texas']), cmap='Reds')
    ax[5].set_title('projected red')

    if pretty:
        ax.axis('off')
        ax.set_yticks([])
        ax.set_xticks([])

    fig.tight_layout()

    plt.savefig(file_plot + '.pdf', dpi=600, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
    print("Saved", file_plot)


def vis_mask(set):
    title = 'mask'
    print("Plotting ", set.file_dat)
    data = np.load(set.file_dat + '.npz')

    fig, ax = plt.subplots(frameon=False)

    selem = morph.disk(5)
    thick_skeleton = morph.dilation(data['skeleton'], selem)
    image = data['mask'] + thick_skeleton
    image = np.ma.masked_where(image == 0, image)

    levels = [0, 1, 2,3]
    colors = ['white', '#8d8dbc' , '#51517c']
    cmap, norm = mpl.colors.from_levels_and_colors(levels, colors)

    ax.imshow(image, cmap=cmap, norm=norm, interpolation='none')

    scaling = set.pixel_scaling
    # correction of scalebar_size to get 1,10,50... um
    # used in pretty mode
    allowed_scalebars = np.append([10**i for i in range(10)],
                                [5* 10**i for i in range(10)] )
    scalebar_size = closest_number(allowed_scalebars,
                                    400 * scaling)
    pixels_scalebar = scalebar_size/scaling
    scalebar = AnchoredSizeBar(ax.transData,
                       pixels_scalebar, str(scalebar_size) + r' $\mu$m',
                       'lower right',
                       pad=0.3,
                       color='black',
                       frameon=False,
                       size_vertical=1)


    ax.add_artist(scalebar)

    ax.axis('off')
    ax.set_yticks([])
    ax.set_xticks([])
    # ax.set_title(title)

    fig.tight_layout()
    plt.savefig(set.file_plot + title + '.pdf', dpi=600, bbox_inches='tight')
    plt.close()


def vis_clean(set, key, cmap, thick=False, cbar=True):
    title = key
    print("Plotting ", set.file_dat)
    data = np.load(set.file_dat + '.npz')

    fig, ax = plt.subplots(frameon=False)


    image = data[key]
    if thick:
        selem = morph.disk(5)
        image = morph.dilation(image, selem)
    image = np.ma.masked_where(image == 0, image)


    im = ax.imshow(image, cmap=cmap)

    scaling = set.pixel_scaling
    # correction of scalebar_size to get 1,10,50... um
    # used in pretty mode
    allowed_scalebars = np.append([10**i for i in range(10)],
                                [5* 10**i for i in range(10)] )
    scalebar_size = closest_number(allowed_scalebars,
                                    400 * scaling)
    pixels_scalebar = scalebar_size/scaling
    scalebar = AnchoredSizeBar(ax.transData,
                       pixels_scalebar, str(scalebar_size) + r' $\mu$m',
                       'lower right',
                       pad=0.3,
                       color='black',
                       frameon=False,
                       size_vertical=1)


    ax.add_artist(scalebar)
    ax.axis('off')

    ax.set_yticks([])
    ax.set_xticks([])
    # ax.set_title(title)

    cbar = fig.colorbar(im)
    fig.tight_layout()

    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('intensity (a.u.)', rotation=270)

    plt.savefig(set.file_plot + title + '.pdf', dpi=600, bbox_inches='tight')
    plt.close()




def plot_samples(data_sets, show, pretty=False):
    for set in data_sets:
        print(set.file_dat)
        print(set.lower_thresh)
        processing_visualization(set.file_dat, set.file_plot, set.analyse_flow,
                                    show=show, pretty=pretty)
        ##

        vis_mask(set)

        vis_clean(set, 'green_clean', 'Greens', thick=False)
        vis_clean(set, 'c_green', 'Greens', thick=True)
        vis_clean(set, 'c_texas', 'Reds', thick=True)





def main():

    set_keyword = os.sys.argv[1].strip()

    start       = int(os.sys.argv[2])
    stop        = int(os.sys.argv[3])
    sample_step = int(os.sys.argv[4])

    plt_range = np.arange(start, stop, sample_step)

    method = 'inter_mean'
    color = 'sep'
    data_sets = [data(set_keyword, i, method, color) for i in plt_range]

    plot_samples(data_sets, False)


if __name__ == "__main__":
    main()
