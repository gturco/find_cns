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
coexp = np.memmap(path + '/coexp.bin', dtype=np.float32, 
                  shape=(len(idxs), len(idxs)), mode='r')

def get_pair(pair):
    pair = pair.strip()
    if not pair: return None, None, None
    if "," in pair: sep = ","
    elif "\t" in pair: sep = "\t"
    else: sep = " "
    pair = [p.strip().upper().replace("_", "") for p in pair.split(sep)]
    assert len(pair) == 2, pair
    return pair, idxs.get(pair[0], -1), idxs.get(pair[1], -1)

REFERENCE = """<p>
<a href="http://nar.oxfordjournals.org/cgi/content/full/37/suppl_1/D987">Reference</a>
<p>
    methods: <a href="http://atted.jp/help/coex_cal.shtml">coexpression</a> and <a href="http://atted.jp/help/mr.shtml">mutual rank</a>
<p>
(the re-scaled value reported here is (1 - MR) / L. where mr is their mutual rank and
L is the number of probes. So a value of 1 means perfect coexpression.
</p>
"""

class form(object):
    def GET(self, data=None):
        web.header('Content-type', 'text/html')
        form = """
        <form method="POST" action=''>
        <p>enter pairs below:</p>
        <p>e.g: AT1G47770,AT2G26550</p>
        <textarea name="pairs" rows="30" cols="20"></textarea>
        <input type="submit" name="submit" />
        </form>
        """
        if data:
            form += "<textarea cols='50' rows='%i'>%s</textarea>" % (max(len(data), 50), "\n".join(data))
        form += REFERENCE
        return form

    def POST(self):
        web.header('Content-type', 'text/plain')
        input_pairs = web.input().pairs.strip().split("\r\n")
        xs = []
        ys = []
        pairs = []
        for pair in input_pairs:
            pair, aidx, bidx = get_pair(pair)
            if pair is None: continue
            xs.append(aidx) 
            ys.append(bidx) 
            pairs.append((pair, -1 in (aidx, bidx)))
        vals = coexp[xs, ys].tolist()
        data = ["#a,b,coexp"]
        for (pair, bad), val in zip(pairs, vals):
            pair = ",".join(pair)
            if bad:
                data.append( "%s,notfound" % pair)
            else:
                data.append( "%s,%.5f" % (pair, val))
        return self.GET(data)


class query(object):
    def GET(self, a, b):
        web.header('Content-type', 'text/plain')
        aidx = idxs.get(a.upper().replace("_", ""))
        bidx = idxs.get(b.upper().replace("_", ""))
        if aidx is None:
            return "%s not found" % aidx
        elif bidx is None:
            return "%s not found" % bidx
        return str(coexp[aidx][bidx])

application = web.application(urls, globals()).wsgifunc()
