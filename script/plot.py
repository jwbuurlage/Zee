#!/usr/bin/python3

# plot.py reads descriptive.mtx, .plt, ... files and plots these` using matplotlib
#
# FIXME: REQUIRES USETEX, PNGDVI, etc.
# TODO: zplot support

import argparse
import os
import yaml

import matplotlib
#matplotlib.use('GTK3Cairo')

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import numpy as np

from math import log, ceil

parser = argparse.ArgumentParser(description=
    "This script reads one or multiple .mtx file(s), and outputs"
    " a spy plot to screen or as a .pdf file")

parser.add_argument('--save', action='store_true',
        help="save to file, dont show on screen")

parser.add_argument('--showfile', action='store_true',
        help="save to file, show file on screen")

parser.add_argument('--filetype', type=str, default='pdf',
        help="filetype used for saving images")

parser.add_argument('--directory', type=str, default='',
        help="the directory in which the matrices are stored. Including"
        " the trailing /")

parser.add_argument('files', type=str, nargs='+',
        help="The file(s) to use as input")

args = parser.parse_args()


def get_color(proc):
    # ith element of has denumerator of 2^ceil(2log(p) + 1)
    #                    numerator is n - whats above, 1 + that
    level = ceil(log(proc + 1, 2))
    denumerator = 2**level
    numerator = 1 + (proc - (denumerator / 2))*2
    color = numerator / denumerator
    color -= 0.5
    if (color < 0.0):
        color += 1.0
    return matplotlib.colors.hsv_to_rgb([color, 1, 1])

def finalize_plt(filename):
    extension_length = len(filename.split('.')[-1])
    if args.save:
        plt.savefig(filename[:-extension_length] + args.filetype)
    elif args.showfile:
        outfilename = filename[:-extension_length] + args.filetype
        plt.savefig(outfilename)
        os.system("xdg-open " + outfilename)
    else:
        plt.show()

###################################################
# SPY MTX

def spy(matrix_file):
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
            finalize_plt(matrix_file)
        elif line.startswith("%%Extended-MatrixMarket distributed-matrix coordinate"):
            while (len(line) == 0 or line[0] == '%'):
                line = fin.readline()
            title = line

            m, n, nz = map(int, fin.readline().split(' '))

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, aspect='equal')

            # force integer labels
            ax.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
            ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))

            # set the appropriate limits
            plt.xlim([0, n])
            plt.ylim([m, 0])

            for _ in range(0, nz):
                i, j, p = map(int, fin.readline().split(' '))
                ax.add_patch(patches.Rectangle((j + marker_offset,
                    i + marker_offset),
                    marker_size, marker_size, color=get_color(p)))

            finalize_plt(matrix_file)
        else:
            print("ERROR: Unrecognized matrix format.")

###################################################
# PLOT YAML

def plot(plot_file):
    f = open(plot_file, 'r')

    # FIXME: check if zee plot file or error
    contents = f.read()
    contents = contents.replace("\\", "\\\\")

    plot_data = yaml.load(contents)

    plt.title(plot_data["title"])
    plt.xlabel(plot_data["xlabel"])
    plt.ylabel(plot_data["ylabel"])
    plt.yscale(plot_data["yscale"])

    for line in plot_data["lines"]:
        line_data = plot_data["lines"][line]["data"]
        plt.plot(line_data)

    finalize_plt(plot_file)

for f in args.files:
    f = args.directory + f
    extension = f.split('.')[-1]
    if extension == "mtx":
        spy(f)
    elif extension == "yaml":
        plot(f)
    else:
        print("can not plot file with extension " + extension)
