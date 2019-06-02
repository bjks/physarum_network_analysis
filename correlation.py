from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from skeleton import *
from analysis.plotting import *


from skimage.feature import register_translation
import os

def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs


def averaged_temporal_correlation(kymo1, kymo2, max_delta, label, path_name):
    squared_kymo1 = crop_aligned_kymo(kymo1)# both kymos have equal shape!
    squared_kymo2 = crop_aligned_kymo(kymo2)

    frame_kymo2 = squared_kymo2[max_delta:-max_delta]
    frame_kymo2 -= np.mean(frame_kymo2)

    t_frame = np.shape(frame_kymo2)[0]
    corr = []
    for i in range(0,2*max_delta):
        temp = squared_kymo1[i:i+t_frame]-np.mean(squared_kymo1[i:i+t_frame])
        corr.append(np.mean((temp) * frame_kymo2))

    t_offset = np.arange(-max_delta,max_delta)
    ######## plotting ########
    plt.ylabel('correlation')
    plt.xlabel('frames offset (T in r(t+T)c(t) )')
    plt.grid(True)

    plt.plot(t_offset, corr)
    plt.savefig(path_name + 'branch' + str(label) + '_corr_av.pdf')
    plt.close()
    return t_offset, corr

##########################################################################################
def temporal_correlation(kymo1, kymo2, max_delta, label, path_name):
    squared_kymo1 = crop_aligned_kymo(kymo1)# both kymos have equal shape!
    squared_kymo2 = crop_aligned_kymo(kymo2)

    frame_kymo2 = squared_kymo2[max_delta:-max_delta]

    correlation = []
    for arr1,arr2 in zip(np.transpose(squared_kymo1), np.transpose(frame_kymo2)):
        arr1-=np.mean(arr1)
        arr2-=np.mean(arr2)
        temp = np.correlate(arr1, arr2, mode = 'valid')
        correlation.append(temp)

    ######## plotting ######
    fig, ax = plt.subplots(1,1, figsize=(3, 6))
    c = ax.imshow(correlation)
    ax.set_ylabel('x')
    ax.set_xlabel('frames offset (T in r(t+T)c(t) )')

    x_range, t_range = np.shape(correlation)[0], np.shape(correlation)[1]
    plt.xticks(np.arange(0,t_range+1, t_range/2), (-t_range/2 ,0, t_range/2))

    plt.colorbar(c)
    plt.grid(True)

    plt.savefig(path_name + 'branch' + str(label) + '_corr_temp.pdf')
    plt.close()
    return np.arange(-max_delta,max_delta+1), correlation

##########################################################################################




######################## PHASE CORR ###########################
def phase_corr(kymo1, kymo2,
               title1, title2, label, path_name,
               align_keyword, alignment, detrending='gauss', show=False, upsample=1):

    kymo1 = align_kymo(kymo1, align_keyword, alignment=alignment)
    kymo2 = align_kymo(kymo2, align_keyword, alignment=alignment)

    squared_kymo1 = np.transpose(crop_aligned_kymo(kymo1))
    squared_kymo2 = np.transpose(crop_aligned_kymo(kymo2))

    if detrending == 'gauss':
        squared_kymo1 = gauss_detrend(squared_kymo1, 40)
        squared_kymo2 = gauss_detrend(squared_kymo2, 40)
    elif detrending == 'mean':
        squared_kymo1 -= np.mean(squared_kymo1)
        squared_kymo2 -= np.mean(squared_kymo2)
    elif detrending == 'radius':
        squared_kymo2/= (squared_kymo1)
        squared_kymo1 = gauss_detrend(squared_kymo1, 30)
        squared_kymo2 = gauss_detrend(squared_kymo2, 30)

    ######## plot kymos ######
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_title(title1)
    ax.set_xlabel("time (frames)")
    ax.set_ylabel("space (pixel)")

    r = ax.imshow(squared_kymo1)
    cbar = plt.colorbar(r)
    cbar.set_label('radius deviations (pixel)', rotation=270)
    plt.savefig(path_name + 'branch' + str(label) + '_' + title1 + '_' +detrending + '.pdf', dpi=600)
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_title(title2)
    ax.set_xlabel("time (frames)")
    ax.set_ylabel("space (pixel)")
    c = ax.imshow(squared_kymo2)
    cbar = plt.colorbar(c)
    cbar.set_label('concentration (a.u.)', rotation=270)

    plt.savefig(path_name + 'branch' + str(label) + '_' + title2 + '_' +detrending + '.pdf', dpi=600)
    plt.close()
    ######

    #### correlation
    squared_kymo1 = ndi.zoom(squared_kymo1, upsample, order=5)
    squared_kymo2 = ndi.zoom(squared_kymo2, upsample, order=5)
    image_product = np.fft.fft2(squared_kymo1) * np.fft.fft2(squared_kymo2).conj()
    cc_image = np.fft.ifft2(image_product)

    unshifted_correlation = cc_image.real
    correlation = np.fft.fftshift(cc_image).real


    # shift, error, diffphase = register_translation(squared_kymo1, squared_kymo2)

    plt.imshow(correlation)

    plt.ylabel('space lag')
    plt.xlabel('time lag')


    x_range, t_range = np.shape(correlation)[0], np.shape(correlation)[1]
    max = np.argmax(unshifted_correlation[0])
    min = np.argmin(unshifted_correlation[0])
    print('Max: ', max/upsample, 'Min: ', min/upsample)
    if min > t_range/2:
        min -= t_range
    if max > t_range/2:
        max -= t_range

    plt.axvline(x=t_range/2 + max, linewidth=0.1, color='k', label='Max: '+str(max/upsample))
    plt.axvline(x=t_range/2 + min, linewidth=0.1, color='r', label='Min: '+str(min/upsample))

    plt.xticks([t_range/2], (0,) )
    plt.yticks([x_range/2], (0,) )

    plt.colorbar()
    plt.grid(linestyle='-', linewidth=0.1)
    plt.legend()

    plt.savefig(path_name + 'branch' + str(label) + '_' + title1 + '_' + title2 +  '_corr_full.pdf', dpi=600)
    if show:
        plt.show()
    plt.close()

def correlate_the2(kymos, quan1, quan2, title1, title2, align_keyword, label,
                   path_name, detrending='gauss', show=False, upsample=1):

    alignment=kymos['alignment']
    kymo1 = kymos[quan1]
    kymo2 = kymos[quan2]

    plot_kymos(kymo1, kymo2, title1, title2, label, path_name, align_keyword, alignment, show=show)
    phase_corr(kymo1, kymo2, title1, title2, label, path_name, align_keyword,
                alignment, detrending=detrending, upsample=upsample, show=show)




def main():
    set_keyword     = os.sys.argv[1]
    color           = os.sys.argv[2]

    align_keyword   = 'reference_point'
    method          = 'inter_mean'


    if color.strip() == 'both':
        colors = ['tg', 'gt']
    else:
        colors = [color]

    labels = range(len(data(set_keyword).seed_positions))
    print(labels[-1])
    for order in colors:
        for label in labels:
            set     = data(set_keyword, method=method, color=order)
            kymos = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')
            path_name = set.file_plot_set + '_branch_' + str(label) + '/'
            if not os.path.exists(path_name):
                os.mkdir(path_name)

            correlate_the2(kymos, 'kymograph_local_radii', 'kymograph_concentration',
                           'radius', 'Ca-concentration',
                           align_keyword, label, path_name, upsample=5, detrending='gauss', show=False)
    #
    # correlate_the2(kymos, 'kymograph_inner', 'kymograph_outer',
    #                'inner-Ca', 'outer-Ca',
    #                align_keyword, label, path_name, 5, show=False)
    #
        # correlate_the2(kymos, 'kymograph_local_radii', 'kymograph_inner',
        #                'radius', 'inner-Ca',
        #                align_keyword, label, path_name, upsample=1, detrending='gauss', show=False)

    # ######### plot kymos ########
    # kymograph_local_radii   = kymos['kymograph_local_radii']
    # kymograph_concentration = kymos['kymograph_concentration']
    # kymograph_inner         = kymos['kymograph_inner']
    # kymograph_outer         = kymos['kymograph_outer']
    #
    # plot_kymos(kymograph_local_radii, kymograph_concentration, 'radii', 'concentration',
    #            label, path_name, align_keyword, alignment)
    #
    # plot_kymos(kymograph_inner, kymograph_outer, 'innner', 'outer',
    #            label, path_name, align_keyword, alignment)
    #
    #
    # ########## correlation stuff
    # radii = align_kymo(kymograph_local_radii, align_keyword, alignment=alignment)
    # conce = align_kymo(kymograph_concentration, align_keyword, alignment=alignment)
    # inner = align_kymo(kymograph_inner, align_keyword, alignment=alignment)
    # outer = align_kymo(kymograph_outer, align_keyword, alignment=alignment)
    #
    # # averaged_temporal_correlation(radii, conce, max_offset, label, path_name)
    # # temporal_correlation(radii, conce, max_offset, label, path_name)
    #
    # phase_corr(radii, conce, 'radii', 'concentration', label, path_name, upsample=30)
    # phase_corr(inner, outer, 'inner', 'outer', label, path_name, upsample=30)



if __name__ == '__main__':
    main()
