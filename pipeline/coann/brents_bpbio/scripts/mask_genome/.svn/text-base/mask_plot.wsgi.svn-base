#!/usr/bin/python
import sys
sys.stdout = sys.stderr
from cgi import parse_qsl
import tables
import os
os.environ[ 'HOME' ] = '/tmp/'
import matplotlib
matplotlib.interactive(0)
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from cStringIO import StringIO

DIR = os.path.abspath(os.path.dirname(__file__))
f = 'copy_count.h5'
h5 = tables.openFile(os.path.join(DIR, f))

def application(env, start_response):
    start_response("200 OK", [("Content-type", "image/png")])

    p = dict(parse_qsl(env['QUERY_STRING']))
    org = p['org'] # e.g. thaliana_v7


    xmin = int(p['xmin'])
    xmax = int(p['xmax'])

    seqid = "c%i" % int(p['seqid'])
    xmin0 = max(xmin, 0)
    vals = h5.getNode('/%s/%s' % (org, seqid))[xmin0:xmax+1]
    ymax = 200

    io = StringIO()

    f = Figure(frameon=False)
    f.canvas = FigureCanvas(f)
    dpi = 128.
    f.set_size_inches(int(p['width'])/dpi, 300/dpi)
    ax = f.add_axes((0,0,1,1), frameon=False, xticks=(), yticks=())
    if xmax > 0:
        ax.set_autoscale_on(0)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(0, ymax)
        ax.plot(range(xmin0, xmax+1), vals)
        ax.plot([xmin0, xmax],[50, 50],'r')

    f.savefig(io, format='png', dpi=dpi)
    io.seek(0)

    return [io.read()]
