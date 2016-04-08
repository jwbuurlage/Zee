#! /usr/bin/python3

from matplotlib import pyplot as pp
from matplotlib import style
from matplotlib2tikz import save as tikz_save

import numpy as np
from math import log, ceil

import numpy as np

fig = pp.figure()
# style.use('ggplot')

D = np.linspace(-0.99, 0.99, 50, True);
I = list(map(lambda x: log((1+x)/(1-x)), D))


pp.xlabel('$x$')
pp.ylabel('$\\log \\left( \\frac{1 + x}{1 - x} \\right)$')

pp.title('behaviour of $w_{\\text{HLP}}$')
pp.grid(True)


pp.plot(D, I, '-', lw=6.0)

tikz_save('../tex/tikz/fn.tikz', figureheight='8cm', figurewidth='12cm');
