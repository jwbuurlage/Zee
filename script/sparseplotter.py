import os
import yaml
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import numpy as np

from math import log, ceil

def get_color(proc):
    # ith element of has denumerator of 2^ceil(2log(p) + 1)
    level = ceil(log(proc + 1, 2))
    denumerator = 2**level
    numerator = 1 + (proc - (denumerator / 2)) * 2
    color = numerator / denumerator
    color -= 0.5
    if (color < 0.0):
        color += 1.0
    return matplotlib.colors.hsv_to_rgb([color, 1, 1])


def finalize_plt(filename):
    extension_length = len(filename.split('.')[-1])
    if args.filetype == "tikz":
        print("INFO: Saved:", filename[:-extension_length] + args.filetype)
        tikz_save(filename[:-extension_length] + args.filetype)
    elif args.save or args.showfile:
        print("INFO: Saved:", filename[:-extension_length] + args.filetype)
        outfilename = filename[:-extension_length] + args.filetype
        plt.savefig(outfilename, bbox_inches='tight')
        if args.showfile:
            os.system(args.viewer + " " + outfilename)
    else:
        print("INFO: Showing")
        plt.show()


def spy(matrix_file, size=3):
    marker_size = 0.8
    marker_offset = 0.5 * (1 - marker_size)

    print("INFO: Opening:", matrix_file)
    with open(matrix_file, 'r') as fin:
        line = fin.readline()
        if line.startswith("%%MatrixMarket matrix array"):
            while (len(line) == 0 or line[0] == '%'):
                line = fin.readline()
            m, n = map(int, line.split(' '))
            data = np.ndarray((m, n))
            for i in range(0, m):
                for j in range(0, n):
                    data[i][j] = float(fin.readline())
            plt.imshow(data)
        elif line.startswith("%%MatrixMarket matrix coordinate"):
            while (len(line) == 0 or line[0] == '%'):
                line = fin.readline()
            title = line

            m, n, nz = map(int, fin.readline().split(' '))

            fig = plt.figure(figsize=(size, size))
            ax = fig.add_subplot(1, 1, 1, aspect='equal')

            # force integer labels
            ax.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
            ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
            plt.xticks([])
            plt.yticks([])
            # plt.axis('off')

            # set the appropriate limits
            plt.xlim([0, n])
            plt.ylim([m, 0])

            for _ in range(0, nz):
                stats = fin.readline().split(' ')
                i, j, p = map(int, stats)
                ax.add_patch(patches.Rectangle((j + marker_offset,
                                                i + marker_offset),
                                               marker_size, marker_size,
                                               color=get_color(p)))
        else:
            print("ERROR: Unrecognized matrix format.")

    plt.show()
