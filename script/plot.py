#!/usr/bin/python3

# plot.py reads descriptive.mtx, .plt, ... files and plots these` using matplotlib
#
# FIXME: REQUIRES USETEX, PNGDVI, etc.
# TODO: zplot support

import argparse
import os
import yaml

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import numpy as np

from sparseplotter import *

#from matplotlib2tikz import save as tikz_save

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

parser.add_argument('--figsize', type=int, default=3,
        help="size in inches of the figure")

parser.add_argument('--viewer', type=str, default='xdg-open',
        help="the viewer for the image")

parser.add_argument('files', type=str, nargs='+',
        help="The file(s) to use as input")

args = parser.parse_args()

###################################################
# PLOT YAML

#def plot(plot_file):
#    f = open(plot_file, 'r')
#
#    # FIXME: check if zee plot file or error
#    contents = f.read()
#    contents = contents.replace("\\", "\\\\")
#
#    plot_data = yaml.load(contents)
#
#    attributes = ["title", "xlabel", "ylabel", "yscale"]
#    for attr in attributes:
#        if attr in plot_data:
#            getattr(plt, attr)(plot_data[attr])
#
#    for line in plot_data["lines"]:
#        line_data = plot_data["lines"][line]["data"]
#        plt.plot(line_data)
#
#    finalize_plt(plot_file)

for f in args.files:
    f = args.directory + f
    extension = f.split('.')[-1]
    if extension == "mtx":
        spy(f)
    elif extension == "yaml":
        plot(f)
    else:
        print("can not plot file with extension " + extension)
