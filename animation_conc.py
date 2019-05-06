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
step = int(os.sys.argv[2])
method = 'inter_mean'

data_sets = [data(set_keyword, i, method) for i in np.arange(data(set_keyword).first, data(set_keyword).last)]

#################################################################


map = []
print('>> Reading...')
first, last = 250, 450
t_arr = np.arange(first,last, step)

for i in t_arr:
    print(' >>', i )
    concentration = np.load(data_sets[i].file_dat + '.npz')['concentration']
    mask = np.where(np.load(data_sets[i].file_dat + '.npz')['green_clean'] !=0, 1, 0)
    skeleton = np.load(data_sets[i].file_dat + '.npz')['skeleton']
    green_clean = np.load(data_sets[i].file_dat + '.npz')['green_clean']
    texas_clean = np.load(data_sets[i].file_dat + '.npz')['texas_clean']


    if keyword == 'raw':
        image = calc_ratio(green_clean, texas_clean)

    else:
        image = tube_radius_at_point(mask, skeleton, concentration)

    image = np.ma.masked_where(image == 0, image)
    map.append(image)

print('>> Create animation...')
fig = plt.figure()
im = plt.imshow(map[0], interpolation="none", cmap='jet', vmin=0.5, vmax=1.2)
plt.colorbar(im)

title = plt.title("")

ani = animation.FuncAnimation(fig, func=update, frames=len(t_arr),
                              repeat=False, interval=200)

print('>> Saving...')
ani.save(data_sets[i].file_plot_set + keyword + '2.mp4', writer="ffmpeg")
