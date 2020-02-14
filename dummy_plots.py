"""

"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
import matplotlib.cm as mc
import matplotlib.colors as mcol

import numpy as np




def main(): ## python3 wall_bf.py <keyword>
    c = np.linspace(0.5,1.5, 100)
    kappa_c = 100

    rate1 = c - 10*(c-1.)**3
    rate2 = c - 50*(c-1.)**3
    rate3 = c - 100*(c-1.)**3




    plt.plot(c, rate1, label = '10' )
    plt.plot(c, rate2, label = '50' )
    plt.plot(c, rate3, label = '100' )

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
