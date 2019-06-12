import numpy as np
import matplotlib.pyplot as plt



def main():

    fig, ax = plt.subplots()

    colors = ['tg', 'gt']
    plot_color = ['blue', 'orange']
    for color, pc in zip(colors, plot_color):
        data_sets_summary = 'time_shifts_data_sets/time_shift_' + color +'.txt'


        set_keywords = np.loadtxt(data_sets_summary, usecols=(0,), dtype='str')
        cols = (6,8)
        extrema_k, extrema_p = np.loadtxt(data_sets_summary, usecols=cols, unpack=True)


        indx = np.argsort(set_keywords)

        set_keywords = set_keywords[indx]
        extrema_k   = extrema_k[indx]

        ax.scatter(range(len(extrema_k)), extrema_k, color=pc, label='order: '+ color)



    ####### plotting #######

    ax.set_ylabel('time_shift')
    ax.set_xticks(range(len(extrema_k)))
    ax.set_xticklabels(set_keywords)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
