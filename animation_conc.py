import matplotlib.animation as animation
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
import os

def update(t):
    im.set_array(map[t])
    title.set_text( 'frame: ' + str(t_arr[t]) )

#################################################################
keyword = 'raw'
title = 'intensity'
step = 1

set_keyword = os.sys.argv[1]
if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
    color = 'gt'
else:
    color = os.sys.argv[3]

step = int(os.sys.argv[2])
method = 'inter_mean'

data_sets = [data(set_keyword, i, method, color=color) for i in np.arange(data(set_keyword).first, data(set_keyword).last)]

#################################################################


map = []
print('>> Reading...')
first, last = data(set_keyword).first , data(set_keyword).last
t_arr = np.arange(first,last, step)
print(t_arr)

for i in range(len(t_arr)-1):
    print(' >>', i )
    print(data_sets[i].file_dat)
    concentration = np.load(data_sets[i].file_dat + '.npz')['concentration']
    mask = np.load(data_sets[i].file_dat + '.npz')['mask']
    skeleton = np.load(data_sets[i].file_dat + '.npz')['skeleton']



    if keyword == 'raw':
        green_clean = np.load(data_sets[i].file_dat + '.npz')['green_clean']
        texas_clean = np.load(data_sets[i].file_dat + '.npz')['texas_clean']
        image = calc_ratio(green_clean, texas_clean)

    else:
        image = tube_radius_at_point(mask, skeleton, concentration)

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
ani.save(data_sets[i].file_plot_set + keyword + '2.mp4', writer="ffmpeg")
