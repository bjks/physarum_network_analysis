import matplotlib.pyplot as plt
import matplotlib.colors as colors
from analysis.network_analysis import *
from analysis.tools import *

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar # scalebar


def plot_image(file_dat, keyword, file_plot, show=False):
    image = np.load(file_dat + '.npz')[keyword]

    fig, ax = plt.subplots()
    image = np.ma.masked_where(image == 0, image)

    c = ax.imshow(image)
    plt.colorbar(c)

    plt.savefig(file_plot + keyword + '.pdf', dpi=600)
    plt.close()


def plot_branch(branch_map, label, path_name, scaling=None, back=None):
    fig, ax = plt.subplots(frameon=False)

    dummy_cmap = colors.ListedColormap(['tab:orange', 'tab:blue'])

    if back!=None:
        ax.imshow(read_file(back), cmap='Greys_r')

    image = -thick_skeleton(branch_map)
    image = np.ma.masked_where(image == 0, image)
    ax.imshow(image, cmap=dummy_cmap)

    if scaling!=None:
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
    plt.savefig(path_name + 'branch' + str(label) + '.pdf', dpi=600, bbox_inches='tight')
    plt.close()
