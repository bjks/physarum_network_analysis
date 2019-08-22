from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *
from ratiometric import *
from skeleton import *
import matplotlib.pyplot as plt


def processing_visualization(file_dat, file_plot, show=False):
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(2,3, figsize=(6, 6), sharex=True, sharey=True)
    ax = axes.ravel()
    green           = data['green_clean']
    green = np.ma.masked_where(green == 0, green)

    texas           = data['texas_clean']
    texas = np.ma.masked_where(texas == 0, texas)

    ax[0].imshow(green, cmap='Greens') #, vmin = min, vmax = max)
    ax[0].set_title('green channel')

    ax[1].imshow(texas, cmap='Reds')
    ax[1].set_title('red channel')

    ax[2].imshow(thick_skeleton(data['skeleton']))
    ax[2].set_title('skeleton')

    ax[3].imshow(thick_skeleton(data['local_radii']))
    ax[3].set_title('mapped radius')

    ax[4].imshow(thick_skeleton(data['c_green']), cmap='Greens')
    ax[4].set_title('projected green')

    ax[5].imshow(thick_skeleton(data['c_texas']), cmap='Reds')
    ax[5].set_title('projected texas')

    fig.tight_layout()
    plt.savefig(file_plot + '.pdf', dpi=600)
    if show:
        plt.show()
    plt.close()



def plot_samples(data_sets, show):
    if not os.path.exists(data_sets[0].path_plots):
        os.mkdir(data_sets[0].path_plots)

    for set in data_sets:
        processing_visualization(set.file_dat, set.file_plot, show=show)


def main():

    set_keyword = os.sys.argv[1].strip()
    sample_step = int(os.sys.argv[4])
    method = 'inter_mean'
    order = 'sep'
    data_sets = [data(set_keyword, i, method, order) for i in np.arange( int(os.sys.argv[2]),int(os.sys.argv[3]), sample_step )]

    plot_samples(data_sets, False)


if __name__ == "__main__":
    main()
