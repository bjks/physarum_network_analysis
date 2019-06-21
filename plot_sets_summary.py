import numpy as np
import matplotlib.pyplot as plt
import os



def main():
    colors           = os.sys.argv[2:]

    fig, ax = plt.subplots()

    means = []

    plot_color = ['blue', 'orange']
    for color, pc in zip(colors, plot_color):
        data_sets_summary = 'time_shifts_data_sets/time_shift_' + color +'.txt'

        set_keywords = np.loadtxt(data_sets_summary, usecols=(0,), dtype='str')

        if os.sys.argv[1]=='phase':
            if color=='green':
                cols = (9,)
            else:
                cols = (8,)
        else:
            if color=='green':
                cols = (7,)
            else:
                cols = (6,)

        extrema_k = np.loadtxt(data_sets_summary, usecols=cols, unpack=True)


        indx = np.argsort(set_keywords)

        set_keywords = set_keywords[indx]
        extrema_k   = extrema_k[indx]
        mean_k = np.mean(extrema_k)
        means.append(mean_k)

        ax.scatter(range(len(extrema_k)), extrema_k, color=pc, label='order: '+ color)
        ax.plot([0, len(extrema_k)], [mean_k, mean_k], color=pc, linestyle='--')


    if len(colors)>1:
        mean = np.mean(means)
        ax.plot([0, len(extrema_k)], [mean, mean], color='purple',
                linestyle='--', label='mean '+str(mean))

    ax.set_ylabel('time shift (s)')
    ax.set_xticks(range(len(extrema_k)))
    ax.set_xticklabels(set_keywords, rotation=30)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
