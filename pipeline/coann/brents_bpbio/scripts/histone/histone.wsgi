#!/usr/bin/python
import sys
sys.stdout = sys.stderr
from cgi import parse_qsl
import os
os.environ[ 'HOME' ] = '/tmp/'
import matplotlib
matplotlib.interactive(0)
matplotlib.rcParams['path.simplify'] = True

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from cStringIO import StringIO
import numpy as np

path = os.path.dirname(__file__)
histones = {}
for ichr in range(1, 6):
    schr = str(ichr)
    f = os.path.join(path, "histone.%s.bin") % schr
    assert os.access(f, os.R_OK)

    histones[schr] = np.memmap(open(os.path.join(path, 'histone.%s.bin') % schr), dtype=np.float32, mode='r')

colors = [(1, 0, 0), (1, 1, 0), (0, 0, 1)]


def application(env, start_response):
    start_response("200 OK", [("Content-type", "image/png")])

    p = dict(parse_qsl(env['QUERY_STRING']))

    xmin = int(p['xmin'])
    xmax = max(int(p['xmax']), 0)
    xmin1 = max(1, xmin)

    seqid = p['seqid']


    hist = np.array(histones[seqid][xmin - 1: xmax])

    height = int(p.get('height', 300))

    f = Figure(frameon=False)
    f.canvas = FigureCanvas(f)
    dpi = 128.
    f.set_size_inches(int(p['width'])/dpi, height/dpi)

    ax = f.add_axes((0, 0, 1, 1), frameon=False, xticks=(), yticks=())
    if xmax != 0:
        ax.set_autoscale_on(0)
        ax.set_xlim(xmin - 1,xmax)

        ax.set_ylim(0.3, 0.9)
        c = 'k' 
        ax.fill_between(
            np.arange(xmin1 - 1, xmax), 
            np.zeros_like(hist),
            hist, facecolor=c, edgecolor='none', alpha=0.4) #, lw=0.8)

    io = StringIO()
    f.savefig(io, format='png', dpi=dpi)
    io.seek(0)
    return [io.read()]
