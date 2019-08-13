import matplotlib.animation as animation
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
import os

def update(t):
    im.set_array(map[t])
    title.set_text( 'frame: ' + str(t_arr[t]) )

#################################################################

set_keyword = os.sys.argv[1]
keyword = os.sys.argv[2]
color = os.sys.argv[3]
step = int(os.sys.argv[4])

method = 'inter_mean'

first, last = data(set_keyword).first, data(set_keyword).last
t_arr = np.arange(first,last, step)

data_sets = [data(set_keyword, i, method, color=color) for i in t_arr]

#################################################################


map = []
print('>> Reading...')
# first, last = data(set_keyword).first , data(set_keyword).last

for set, t in zip(data_sets, t_arr):
    print(' >>', t )
    print(set.file_dat)
    local_radii   = np.load(set.file_dat + '.npz')['local_radii']
    mask = np.load(set.file_dat + '.npz')['mask']
    skeleton = np.load(set.file_dat + '.npz')['skeleton']

    if keyword == 'raw':
        green_clean = np.load(set.file_dat + '.npz')['green_clean']
        texas_clean = np.load(set.file_dat + '.npz')['texas_clean']
        rd = np.load(set.file_dat + '.npz')['rel_dist']

        image = calc_ratio(green_clean, texas_clean)

    elif keyword == 'raw_net':
        image = np.load(set.file_dat + '.npz')['network_clean']

    elif keyword == 'raw_norm':
        rd = np.load(set.file_dat + '.npz')['rel_dist']
        r  = tube_radius_at_point(mask, skeleton, local_radii)
        height = np.where(rd<1, r * np.sqrt(1 - rd**2), 0)
        if color == 'green':
            network = np.load(set.file_dat + '.npz')['network_clean']
        else:
            network = np.load(set.file_dat + '.npz')['green_clean']

        image = np.where(height>0, calc_ratio(network,height), 0)


    elif keyword == 'concentration':
        concentration = np.load(set.file_dat + '.npz')['concentration']

        image = tube_radius_at_point(mask, skeleton, concentration)

    elif keyword == 'norm_conc':
        concentration = np.load(set.file_dat + '.npz')['concentration']

        norm_conc = calc_ratio(concentration, local_radii)
        image = tube_radius_at_point(mask, skeleton, norm_conc)

    elif keyword == 'skeleton':
        image = thick_skeleton(skeleton)

    image = np.where(image == 0, np.nan, image)
    map.append(image)

print('>> Create animation...')
fig = plt.figure()
std_im, mean_im = np.nanstd(map), np.nanmean(map)
print(mean_im)
min, max = mean_im - std_im, mean_im + std_im
im = plt.imshow(map[0], interpolation="none", cmap='Spectral_r', vmin = min, vmax = max)
plt.colorbar(im)

title = plt.title("")

ani = animation.FuncAnimation(fig, func=update, frames=len(t_arr),
                              repeat=False, interval=200)

print('>> Saving...')
if not os.path.exists(data_sets[0].path_plots):
    os.mkdir(data_sets[0].path_plots)
ani.save(data_sets[0].file_plot_set + keyword + '.mp4', writer="ffmpeg")
