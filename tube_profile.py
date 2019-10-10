import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.data_sets import *
from scipy import stats
import os


def circle_height(x):
    return np.sqrt(1. - x**2)

def collect_data(data_sets):
    r   = np.array([])
    d   = np.array([])
    rd  = np.array([])
    g   = np.array([])
    t   = np.array([])

    for set in data_sets:
        print(set.file_dat)

        green_clean     = np.load(set.file_dat + '.npz')['green_clean']
        texas_clean     = np.load(set.file_dat + '.npz')['texas_clean']
        relative_dist   = np.load(set.file_dat + '.npz')['rel_dist']
        radius          = np.load(set.file_dat + '.npz')['radii_map']
        mask            = np.load(set.file_dat + '.npz')['mask']

        dist_mask   = np.where((relative_dist<=1)*(mask!=0), 1, 0)

        r  =    np.append(r,         radius[np.where(dist_mask)])
        rd =    np.append(rd, relative_dist[np.where(dist_mask)])
        g  =    np.append(g,    green_clean[np.where(dist_mask)])
        t  =    np.append(t,    texas_clean[np.where(dist_mask)])

    return r, rd, g, t, calc_ratio(g,t)

def bin_radii(r, rd, g, t, min, max):
    within_interval = np.where( (r>min)*(r<max) )
    return rd[within_interval], g[within_interval], t[within_interval]



def bin_data(r, rd, c, r_bins, rd_bins):

    profile, r_bin_edges, rd_bin_edges,_= stats.binned_statistic_2d(r, rd, c, statistic='mean', bins=[r_bins, rd_bins])
    profile_std, _, _, _                = stats.binned_statistic_2d(r, rd, c, statistic=np.std, bins=[r_bins, rd_bins])

    norm = np.transpose(profile)[0]
    profile     = [p/N for p,N in zip(profile, norm)]
    profile_std = [p/N for p,N in zip(profile_std, norm)]

    return r_bin_edges, rd_bin_edges, profile, profile_std


def main():
    set_keyword = os.sys.argv[1]
    step        = int(os.sys.argv[2])
    color       = 'sep'
    method      = 'inter_mean'

    data_sets = [data(set_keyword, i, color=color, method=method) for i in np.arange(data(set_keyword).first, data(set_keyword).last, step)]
    print(np.arange(data(set_keyword).first, data(set_keyword).last, step))

    radius, rel_dist, green, texas, ratio = collect_data(data_sets)

    min_r = 0.3*np.max(radius)
    max_r = np.max(radius)
    r_bins  = np.around(np.linspace(min_r, max_r, int((max_r-min_r)/2))).astype(int)
    rd_bins = np.arange(0,1.02, 0.02)

    _,_, g, g_std = bin_data(radius, rel_dist, green, r_bins, rd_bins)
    _,_, t, t_std = bin_data(radius, rel_dist, texas, r_bins, rd_bins)
    _,_, r, r_std = bin_data(radius, rel_dist, ratio, r_bins, rd_bins)

    c = plt.pcolormesh(r)
    r_range = r_bins
    plt.yticks(np.arange(len(r_range))[::5], r_range[::5])
    plt.xticks(np.linspace(0, len(rd_bins)-1, 6), [0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.ylabel('radius (pixels)')
    plt.xlabel('relative distance r/R')
    plt.title('ratio')
    plt.colorbar(c)
    filename = data(set_keyword).file_plot_tube_profile + '_grid_tube_profile.pdf'
    plt.savefig(filename, dpi=200)

    plt.show()
    plt.close()

    bin = rd_bins[:-1] + np.diff(rd_bins)/2

    for j in range(len(r_bins) -1):

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(bin, g[j], color='green', label='green channel')
        ax1.fill_between(bin, g[j] - g_std[j], g[j] + g_std[j], alpha=0.2, color='green')
        ax1.plot(bin, t[j], color='red', label='red channel')
        ax1.fill_between(bin, t[j] - t_std[j], t[j] + t_std[j], alpha=0.2, color='red')
        ax2.plot(bin, r[j], color='blue', label='ratio')
        ax2.fill_between(bin, r[j] - r_std[j], r[j] + r_std[j], alpha=0.2, color='blue')

        ax1.plot(rd_bins, circle_height(rd_bins), ':' , color='black')

        ax1.legend()
        # ax1.legend(frameon=False)
        ax1.set_xlabel('relative distance r/R')
        ax1.set_ylabel('normalized profiles')
        ax2.set_ylabel('ratio')
        plt.title(str(r_bins[j]) + '<R<' + str(r_bins[j+1]) + ' pixels')
        filename = data(set_keyword).file_plot_tube_profile + '_' + str(r_bins[j]) + '-' + str(r_bins[j+1]) + '_tube_profile.pdf'
        plt.savefig(filename, dpi=200)
        # plt.show()
        plt.close()

if __name__ == '__main__':
    main()
