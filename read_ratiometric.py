from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *
from ratiometric import *
from skeleton import *
import matplotlib.pyplot as plt


def processing_visualization(file_dat, file_plot, analyse_flow,show=False):
    print("Plooting ", file_dat)
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
        ax[3].imshow(thick_skeleton(data['skeleton']), cmap='Blues_r')



    ax[4].imshow(thick_skeleton(data['c_green']), cmap='Greens')
    ax[4].set_title('projected green')

    ax[5].imshow(thick_skeleton(data['c_texas']), cmap='Reds')
    ax[5].set_title('projected red')



    fig.tight_layout()
    plt.savefig(file_plot + '.pdf', dpi=600, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
    print("Saved", file_plot)




def plot_samples(data_sets, show):
    for set in data_sets:
        processing_visualization(set.file_dat, set.file_plot, set.analyse_flow,
                                    show=show)


def main():

    set_keyword = os.sys.argv[1].strip()

    start       = int(os.sys.argv[2])
    stop        = int(os.sys.argv[3])
    sample_step = int(os.sys.argv[4])

    plt_range = np.arange(start, stop, sample_step)

    method = 'inter_mean'
    order = 'sep'
    data_sets = [data(set_keyword, i, method, order) for i in plt_range]

    plot_samples(data_sets, False)


if __name__ == "__main__":
    main()
