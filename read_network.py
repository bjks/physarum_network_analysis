from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *
from ratiometric import *
from skeleton import *
import matplotlib.pyplot as plt

def processing_visualization_network(file_dat, file_plot, show=False):
    print("Plooting ", file_dat)
    data = np.load(file_dat + '.npz')

    fig, axes = plt.subplots(1,3, sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(data['network_clean'], cmap='Greys_r')
    ax[0].set_title('network')

    radius = np.ma.masked_where(data['radii_map'] == 0, data['radii_map'])
    ax[1].imshow(radius)
    ax[1].set_title('mapped radius')

    ax[2].imshow(thick_skeleton(data['concentration']))
    ax[2].set_title('intensity')

    fig.tight_layout()
    plt.savefig(file_plot + '.pdf', dpi=600, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
    print("Saved", file_plot)



def plot_samples_network(data_sets, show=False):
    for set in data_sets:
        processing_visualization_network(set.file_dat, set.file_plot, show=show)




def main():

    set_keyword = os.sys.argv[1]
    color       = str(os.sys.argv[2])
    start       = int(os.sys.argv[3])
    stop        = int(os.sys.argv[4])
    sample_step = int(os.sys.argv[5])

    method = 'inter_mean'
    plt_range = np.arange(start, stop, sample_step)


    data_sets = [data(set_keyword, i, method, color) for i in plt_range]

    plot_samples_network(data_sets)


if __name__ == "__main__":
    main()
