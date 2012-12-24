#!/usr/bin/python
import web
import numpy as np
import os

urls = (
    '/([\w_\d]+)/([\w_\d]+)/?', 'query',
    '/', 'form'
)

path = os.path.dirname(__file__)

idxs = dict((accn.rstrip(), i) for i, accn in enumerate(open(path + "/names.order")))
aranet = np.memmap(path + '/aranet.bin', dtype=np.float32, 
                  shape=(len(idxs), len(idxs)), mode='r')

def get_pair(pair):
    pair = pair.strip()
    if not pair: return None, None, None
    if "," in pair: sep = ","
    elif "\t" in pair: sep = "\t"
    else: sep = " "
    pair = [p.strip().upper().replace("_", "") for p in pair.split(sep) if p.strip()]
    assert len(pair) == 2, pair
    return pair, idxs.get(pair[0], -1), idxs.get(pair[1], -1)

DEFAULT="""\
AT1G79010,ATMG00270
AT1G16700,ATMG00270
AT2G20580,AT4G19006
AT2G02570,AT3G55220
AT2G02570,AT3G55200
ATMG00900,ATMG00960
AT2G07771,ATMG00960
AT2G07681,ATMG00960
AT4G19006,AT5G64760
"""

REFERENCE = """
<a href="http://www.functionalnet.org/aranet/about.html">Aranet site</a>
<p>
<a href="http://www.nature.com/nbt/journal/v28/n2/abs/nbt.1603.html">Publication</a>
"""

class form(object):
    def GET(self, form_data=DEFAULT, data=None):
        web.header('Content-type', 'text/html')
        form = """
        <form method="POST" action=''>
        <p>enter pairs below:</p>
        <textarea name="pairs" rows="30" cols="20">%s</textarea>
        <input type="submit" name="submit" />
        </form>
        """ % form_data
        if data:
            form += "<textarea cols='50' rows='%i'>%s</textarea>" % (max(len(data), 50), "\n".join(data))
        return form + REFERENCE

    def POST(self):
        web.header('Content-type', 'text/plain')
        pairs_txt  = web.input().pairs.strip()
        input_pairs = pairs_txt.split("\r\n")
        xs = []
        ys = []
        pairs = []
        for pair in input_pairs:
            pair, aidx, bidx = get_pair(pair)
            if pair is None: continue
            xs.append(aidx) 
            ys.append(bidx) 
            pairs.append((pair, -1 in (aidx, bidx)))
        vals = aranet[xs, ys].tolist()
        data = ["#a,b,aranet"]
        for (pair, bad), val in zip(pairs, vals):
            pair = ",".join(pair)
            if bad:
                data.append( "%s,accn-not-found" % pair)
            elif val == 0:
                data.append( "%s,nodata-for-pair" % (pair, ))
            else:
                data.append( "%s,%.5f" % (pair, val))
        return self.GET(form_data=pairs_txt, data=data)


class query(object):
    def GET(self, a, b):
        web.header('Content-type', 'text/plain')
        aidx = idxs.get(a.upper().replace("_", ""))
        bidx = idxs.get(b.upper().replace("_", ""))
        if aidx is None:
            return "%s not found" % aidx
        elif bidx is None:
            return "%s not found" % bidx
        return str(aranet[aidx][bidx])

application = web.application(urls, globals()).wsgifunc()
