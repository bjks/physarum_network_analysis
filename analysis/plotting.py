import matplotlib.pyplot as plt
import matplotlib.colors as colors
from analysis.tools import *
from analysis.network_analysis import *

def plot_image(file_dat, keyword, file_plot, show=False):
    image = np.load(file_dat + '.npz')[keyword]

    fig, ax = plt.subplots()
    image = np.ma.masked_where(image == 0, image)

    c = ax.imshow(image)
    plt.colorbar(c)

    plt.savefig(file_plot + keyword + '.pdf', dpi=600)
    plt.close()


def plot_branch(branch_map, label, path_name):
    dummy_cmap = colors.ListedColormap(['tab:orange', 'tab:blue'])

    plt.imshow(-thick_skeleton(branch_map), cmap=dummy_cmap)
    plt.savefig(path_name + 'branch' + str(label) + '.pdf')
    plt.close()
