import matplotlib.animation as animation
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
import os

def update(t):
    im.set_array(map[t])
    title.set_text(keyword + ': frame ' + str(t_arr[t]) )

#################################################################
keyword = 'disk_mean'
step = 1

set_keyword = os.sys.argv[1]
method = 'disk_mean'

data_sets = [data(set_keyword, i, method) for i in range(data(set_keyword).first, data(set_keyword).last)]
seed_position = data(set_keyword).seed_positions[int(os.sys.argv[2])]

#################################################################


map = []
print('>> Reading...')
first, last = 1, 100
t_arr = np.arange(first,last, 1)

for i in t_arr:
    print(' >>', i )
    concentration = np.load(data_sets[i].file_dat + '.npz')['concentration']
    mask = np.where(np.load(data_sets[i].file_dat + '.npz')['green_clean'] !=0, 1, 0)
    skeleton = np.load(data_sets[i].file_dat + '.npz')['skeleton']
    green_clean = np.load(data_sets[i].file_dat + '.npz')['green_clean']
    texas_clean = np.load(data_sets[i].file_dat + '.npz')['texas_clean']

    seed_position   = closest_point_in_skel(seed_position, skeleton)
    branches        = morph.label(mask, connectivity=2)
    label           = branches[seed_position[0], seed_position[1]]
    mask            = np.where(branches==label, 1, 0)

    if keyword == 'raw':
        green_clean *= mask
        texas_clean *= mask
        image = calc_ratio(green_clean, texas_clean)

    elif keyword == 'disk_mean':
        image = tube_radius_at_point(mask, skeleton, concentration)


    map.append(image)

print('>> Create animation...')
fig = plt.figure()
im = plt.imshow(map[0], interpolation="none", cmap='tab10')
plt.colorbar(im)

title = plt.title("")

ani = animation.FuncAnimation(fig, func=update, frames=len(t_arr),
                              repeat=False, interval=200)

print('>> Saving...')
ani.save(data_sets[i].file_plot_set + keyword + '2.mp4', writer="ffmpeg")
