from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.data_setsBF import *

from analysis.plotting import *
from ratiometric import *
from skeleton import *
import matplotlib.pyplot as plt

def processing_visualization_network(file_dat, file_plot, show=False):
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(1,3, figsize=(6, 6), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(data['network_clean'])
    ax[0].set_title('network')

    ax[1].imshow(thick_skeleton(data['local_radii']))
    ax[1].set_title('Local radii')

    ax[2].imshow(thick_skeleton(data['concentration']))
    ax[2].set_title('Concentration')

    fig.tight_layout()
    plt.savefig(file_plot + '.pdf')
    if show:
        plt.show()
    plt.close()


def plot_samples_network(data_sets, show=True):
    if not os.path.exists(data_sets[0].path_plots):
        os.mkdir(data_sets[0].path_plots)

    for set in data_sets:
        processing_visualization_network(set.file_dat, set.file_plot, show=show)




def main():

    set_keyword = os.sys.argv[1]
    color       =  str(os.sys.argv[2])

    start       = int(os.sys.argv[3])
    stop        = int(os.sys.argv[4])
    sample_step = int(os.sys.argv[5])


    method = 'inter_mean'

    data_sets = [data(set_keyword, i, method, color=color) for i in np.arange(start, stop, sample_step )]

    # data_sets = [data(set_keyword, i, method, order) for i in np.arange(2,3)]

    plot_samples_network(data_sets)


if __name__ == "__main__":
    main()
