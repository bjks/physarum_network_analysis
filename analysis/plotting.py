import matplotlib.pyplot as plt
from  analysis.tools import *
from analysis.network_analysis import *


def processing_visualization(file_dat, file_plot, show=False):
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(2,3, figsize=(6, 6), sharex=True, sharey=True)
    ax = axes.ravel()
    green           = data['green_clean']
    green = np.ma.masked_where(green == 0, green)

    # green_spotless  = data['green_spotless']
    # green = np.ma.masked_where(green == 0, green)
    #
    texas           = data['texas_clean']
    texas = np.ma.masked_where(texas == 0, texas)

    ratio           = data['ratio']
    ratio = np.ma.masked_where(ratio == 0, ratio)


    # max = np.amax([green, green_spotless])
    # min = np.amin([green, green_spotless])


    ax[0].imshow(green, cmap='Greens') #, vmin = min, vmax = max)
    ax[0].set_title('green channel')

    ax[1].imshow(texas, cmap='Reds')
    ax[1].set_title('red channel')

    ax[2].imshow(ratio)
    ax[2].set_title('green/red')

    # ax[1].imshow(green_spotless, cmap='Greens', vmin = min, vmax = max)
    # ax[1].set_title('green without spots')
    #
    # ax[2].imshow(texas, cmap='Reds')
    # ax[2].set_title('red channel')
    #
    # ax[3].imshow(ratio, cmap='Purples')
    # ax[3].set_title('green/red')

    ax[3].imshow(thick_skeleton(data['skeleton']))
    ax[3].set_title('skeleton')

    ax[4].imshow(thick_skeleton(data['local_radii']))
    ax[4].set_title('tube radius')

    ax[5].imshow(thick_skeleton(data['concentration']))
    ax[5].set_title('projected intensity')

    fig.tight_layout()
    plt.savefig(file_plot + '.pdf', dpi=600)
    if show:
        plt.show()
    plt.close()


def processing_visualization_color(file_dat, file_plot, color, show=False):
    file_dat += color
    file_plot +=color
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(1,3, figsize=(6, 6), sharex=True, sharey=True)
    ax = axes.ravel()

    if color == 'green':
        ax[0].imshow(data['green_clean'])
        ax[0].set_title('Green')
    else :
        ax[0].imshow(data['texas_clean'])
        ax[0].set_title('Texas')

    ax[1].imshow(thick_skeleton(data['local_radii']))
    ax[1].set_title('Local radii')

    ax[2].imshow(thick_skeleton(data['concentration']))
    ax[2].set_title('Concentration')

    fig.tight_layout()
    plt.savefig(file_plot + '.pdf')
    if show:
        plt.show()
    plt.close()


def plot_branch(branch_map, label, path_name):
    plt.imshow( thick_skeleton(branch_map))
    plt.savefig(path_name + 'branch' + str(label) + '.pdf')
    plt.close()


def plot_kymos(kymo1, kymo2,
               title1, title2,
               label, path_name, mode, alignment, show=False):

    fig, axes = plt.subplots(1,2, figsize=(10, 4), sharex=True, sharey=True)
    ax = axes.ravel()

    squared_1 = align_kymo(kymo1, mode, alignment)
    r = ax[0].imshow(np.transpose(squared_1))
    ax[0].set_title(title1)
    ax[0].set_ylabel('pixel')
    ax[0].set_xlabel('frames')
    fig.colorbar(r, ax=ax[0])

    squared_2 = align_kymo(kymo2, mode, alignment)
    c = ax[1].imshow(np.transpose(squared_2))
    ax[1].set_title(title2)
    ax[1].set_ylabel('pixel')
    ax[1].set_xlabel('frames')
    fig.colorbar(c, ax=ax[1])

    fig.tight_layout()
    plt.savefig(path_name + 'branch' + str(label) + '_' + title1 + '_' + title2 + '_kymographs.pdf', dpi=600)
    if show:
        plt.show()
    plt.close()
