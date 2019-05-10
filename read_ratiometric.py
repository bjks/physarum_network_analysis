from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.plotting import *
from ratiometric import *
from skeleton import *
import matplotlib.pyplot as plt


def plot_samples(data_sets, show=True):
    if not os.path.exists(data_sets[0].path_plots):
        os.mkdir(data_sets[0].path_plots)

    for set in data_sets:
        processing_visualization(set.file_dat, set.file_plot, show=show)




def main():

    set_keyword = os.sys.argv[1]
    sample_step = int(os.sys.argv[2])
    method = 'inter_mean'
    order = 'gt'
    data_sets = [data(set_keyword, i, method, order) for i in np.arange( int(os.sys.argv[2]),int(os.sys.argv[3]), sample_step )]
    # data_sets = [data(set_keyword, i, method, order) for i in np.arange(2,3)]

    # def plot_image(file_dat, keyword, file_plot, show=False):

    # plot_image(data_sets[0].file_dat, 'ratio', data_sets[0].file_plot)
    plot_samples(data_sets, True)


if __name__ == "__main__":
    main()
