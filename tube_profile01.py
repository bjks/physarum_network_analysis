import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.data_sets import *
from scipy import stats
import os
#
# def calculate_profile(dye, mask, skel, local_radii, step_dist, step_tube):
#     distance = ndi.distance_transform_edt(np.invert(skel.astype(bool)),
#                                       return_distances=True, return_indices=False)
#
#     distance *= mask
#     radii = tube_radius_at_point(mask, skel, local_radii)
#     r_max = np.max(radii)
#
#     relative_dist = np.true_divide(distance, radii, out=np.zeros_like(radii), where=radii!=0)
#
#     bins_dist =     np.arange(0, 1,     step_dist)
#     bins_tube =     np.arange(0, 100,   step_tube)
#
#     profiles = np.zeros( (np.shape(bins_tube)[0], np.shape(bins_dist)[0]) )
#
#     for profile, tube in zip(profiles, bins_tube):
#         for j in range(0, len(bins_dist)):
#             con=(radii > tube)*(radii < tube+step_tube)*(relative_dist > bins_dist[j])*(relative_dist < bins_dist[j]+step_dist)
#             sum     = np.sum(np.where(con, dye, 0))
#             count   = np.sum(np.where(con, 1  , 0))
#             if count != 0:
#                 profile[j] = sum/count
#
#             else:
#                 profile[j] = 0
#     return profiles
#
#
# def concentration_profiles(set, step_dist, step_tube):
#     green = read_file( set.file_raw1 )
#     texas = read_file( set.file_raw2 )
#     mask = create_mask(texas, set.sigma, set.threshold, set.halo_sig)
#
#     green_clean = green * mask
#     texas_clean = texas * mask
#
#
#     ratio                          = calc_ratio(green_clean, texas_clean)
#     skeleton                       = extract_skeleton(mask)
#     local_radii                    = extract_radii(mask, skeleton)
#
#     profiles = calculate_profile(ratio, mask, skeleton, local_radii, step_dist, step_tube)
#     return profiles
#
#
# def plot_profiles(set_keyword, start, stop):
#     # start, stop = int(os.sys.argv[2]), int(os.sys.argv[3])
#     data_sets = [data(set_keyword, i) for i in range(start, stop)]
#
#     step_dist = 0.05
#     step_tube = 10
#
#     profiles = []
#     i = 0
#     for set in data_sets:
#         print('\r>>', i, end='')
#         i+=1
#         profiles.append(concentration_profiles(set, step_dist, step_tube))
#
#
#     mean_profiles = np.mean(profiles, axis = 0)
#     bins_dist =     np.arange(0, 1,     step_dist)
#     bins_tube =     np.arange(0, 100,   step_tube)
#
#     skip = 2
#     for profile, bin_t in zip(mean_profiles[skip:], bins_tube[skip:]):
#         plt.plot(bins_dist, profile, label = str(bin_t) + ' < R < ' + str(bin_t+step_tube))
#
#     plt.xlabel('realtive distance r/R')
#     plt.ylabel('Ca concentration (a.u.)')
#
#
#     plt.legend()
#     filename = set.file_plot_set + '_tube_profile.pdf'
#     plt.savefig(filename, dpi=600)
#     plt.show()

def circle_height(x):
    return np.sqrt(1 - x**2)

def collect_data(data_sets):
    r   = np.array([])
    d   = np.array([])
    rd  = np.array([])
    g   = np.array([])
    t   = np.array([])

    for set in data_sets:
        green = read_file( set.file_raw1 )
        texas = read_file( set.file_raw2 )
        mask = create_mask(texas, set.sigma, set.threshold, set.halo_sig)

        mask = extract_nerwork(mask)

        green_clean = green * mask
        texas_clean = texas * mask

        skeleton                       = extract_skeleton(mask)
        local_radii                    = extract_radii(mask, skeleton)

        relative_dist = relative_distance(skeleton, mask, local_radii)
        radius = tube_radius_at_point(mask, skeleton, local_radii)

        dist_mask   = np.where((relative_dist<=1)*(mask!=0), 1, 0)

        r  =    np.append(r,         radius[np.where(dist_mask)])
        rd =    np.append(rd, relative_dist[np.where(dist_mask)])
        g  =    np.append(g,    green_clean[np.where(dist_mask)])
        t  =    np.append(t,    texas_clean[np.where(dist_mask)])

    return r, rd, g, t, g/t

def bin_radii(r, rd, g, t, min, max):
    within_interval = np.where( (r>min)*(r<max) )
    return rd[within_interval], g[within_interval], t[within_interval]



def bin_data(rd, c, bins):

    profile, bin_edges, binnumber   = stats.binned_statistic(rd, c, statistic='mean', bins=bins)
    profile_std, _, _               = stats.binned_statistic(rd, c, statistic=np.std, bins=bins)

    norm = profile[0]
    profile /= norm
    profile_std /= norm

    return bin_edges, profile, profile_std


def main():
    set_keyword = os.sys.argv[1]
    step        = int(os.sys.argv[2])
    no_bins     = 50
    data_sets = [data(set_keyword, i) for i in np.arange(data(set_keyword).first, data(set_keyword).last, step)]

    r, rd, g, t, ratio = collect_data(data_sets)
    colected_profiles = []

    step = 5
    Rs = np.arange(20,50, step)
    for R in Rs:
        print( R, R+step)
        rd_bin, g_bin, t_bin = bin_radii(r, rd, g,t, R, R+step)

        bin_edges, g_profile, g_std = bin_data(rd_bin, g_bin, no_bins)

        bin = bin_edges[:-1] + np.diff(bin_edges)/2

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(bin, g_profile, color='green', label='green channel')
        ax1.fill_between(bin, g_profile - g_std,
                         g_profile + g_std, alpha=0.2, color='green')

        bin_edges, t_profile, t_std = bin_data(rd_bin, t_bin, no_bins)
        ax1.plot(bin, t_profile, color='red', label='red channel')
        ax1.fill_between(bin, t_profile - t_std,
                         t_profile + t_std , alpha=0.2, color='red')

        bin_edges, r_profile, r_std = bin_data(rd_bin, g_bin/t_bin, no_bins)
        ax2.plot(bin, r_profile, color='blue', label='ratio')
        ax2.fill_between(bin, r_profile - r_std,
                         r_profile + r_std , alpha=0.2, color='blue')
        colected_profiles.append(r_profile)


        ax1.plot(bin_edges, circle_height(bin_edges), ':' , color='black')

        ax1.legend()
        # ax1.legend(frameon=False)
        ax1.set_xlabel('realtive distance r/R')
        ax1.set_ylabel('normalized profiles')
        ax2.set_ylabel('ratio')
        plt.title(str(R) + '<R<' + str(R+step) + ' pixels')
        filename = data(set_keyword).file_plot_set + '_' + str(R) + '-' + str(R+step) + '_tube_profile.pdf'
        plt.savefig(filename, dpi=600)
        plt.close()
##################

    plt.imshow(colected_profiles)
    plt.yticks([0, len(Rs)-1], (0,Rs[-1]) )
    plt.xticks([0, bin_edges[-1]], (0,1) )
    plt.xlabel('average radius (pixels)')
    plt.ylabel('realtive distance r/R')
    plt.colorbar()

    filename = data(set_keyword).file_plot_set + '_coll_' + '_tube_profile.pdf'
    plt.savefig(filename, dpi=600)

    plt.close()


##################
    # fig, axes = plt.subplots(1,3, figsize=(10, 4), sharex=True, sharey=True)
    # ax = axes.ravel()
    #
    # a = 0.1
    # s = 0.5
    # ax[0].scatter(rd, t, alpha=a, s=s, c=r)
    # ax[0].set_title('Texas')
    # ax[0].set_ylabel('intensity')
    # ax[0].set_xlabel('relative distance r/R')
    #
    # ax[1].scatter(rd, g, alpha=a, s=s, c=r)
    # ax[1].set_title('Green')
    # ax[1].set_ylabel('intensity')
    # ax[1].set_xlabel('relative distance r/R')
    #
    # ax[2].scatter(rd, ratio, alpha=a, s=s, c=r)
    # ax[2].set_title('ratio')
    # ax[2].set_ylabel('intensity')
    # ax[2].set_xlabel('relative distance r/R')
    #
    # filename = data(set_keyword).file_plot_set + '_scatter_' + '_tube_profile.pdf'
    # plt.savefig(filename)
    # plt.show()
    # plt.close()
    #

if __name__ == '__main__':
    main()
