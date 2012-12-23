import os
from subprocess import Popen
import sys

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(os.path.basename(__file__))

def get_blast_file(qfasta, sfasta, out_dir=None):
    """
    >>> get_blast_file("/tmp/a.fasta", "b.fasta")
    'a_vs_b.blast'
    >>> get_blast_file("/tmp/a.fasta", "b.fasta", out_dir=True)
    '/tmp/a_vs_b.blast'
    """
    q = os.path.basename(qfasta)
    s = os.path.basename(sfasta)
    blast_file = q[:q.rfind(".")] + "_vs_" + s[:s.rfind(".")] + ".blast"
    if not out_dir: return blast_file
    if out_dir is True:
        d = os.path.dirname(qfasta)
        return os.path.abspath(os.path.join(d, blast_file))
    return os.path.join(out_dir, blast_file)

def is_current_file(a, b):
    """
    >>> is_current_file(__file__, 'a.txt')
    False
    >>> is_current_file(__file__, __file__)
    False
    """
    if not (os.path.exists(a) and os.path.exists(b)): return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm

def is_same_blast_params(blast_file, cmd):
    """ when using the blast() below, a .cmd file is
    written, this function checks that file to see if the
    current command is the same as that. if so, and the fasta
    files are up to date, the blast is not done as it's up to date"""
    params_file = blast_file + ".cmd"
    if not os.path.exists(params_file): return False
    return cmd.strip() == open(params_file).read().strip()

def sh(cmd, blast_log=None):
    """run a commmand in the shell
    # waiting for bugfix in nose.
    #>>> sh("echo 'hi'")
    'hi'
    """
    if not blast_log is None:
        cmd += " 2>%s" % blast_log
    log.debug(cmd)
    proc = Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
    r = proc.communicate()
    return r

def add_dash(params):
    """
    >>> add_dash({'p': 'F', 'save': 'T'})
    {'-p': 'F', '--save': 'T'}
    """
    d = {}
    for k, v in params.items():
        if k.startswith("-"):
            d[k] = v
        elif len(k) > 1:
            d["--" + k] = v
        else:
            d["-" + k] = v
    return d


def is_protein_db(blast_cfg):
    """
    >>> is_protein_db({'p': 'blastx'})
    True
    """
    return blast_cfg["p"] in ("blastx", "blastp")

def rm(f):
    if os.path.exists(f): os.unlink(f)
   

def blast(_blast_cfg, full_name=False, blast_log=None):
    """
    >>> blast({'i': 'tests/a.fasta', 'd': 'tests/a.fasta'})

    """
    blast_cfg = _blast_cfg.copy()
    check_args(blast_cfg)
    q_fasta = blast_cfg["i"]
    s_fasta = blast_cfg["d"]
    blastall = blast_cfg["b"]
    blast_cfg.pop('b')

    format_db = blastall.replace("blastall", "formatdb")

    protein = "T" if is_protein_db(blast_cfg) else "F"
    cmd = "%(format_db)s -i %(s_fasta)s -p %(protein)s" % locals()

    ext = ".pin" if protein == "T" else ".nin"
    assert os.path.exists(s_fasta), "%s does not exist!" % s_fasta
    if not is_current_file(s_fasta + ext, s_fasta):
        try:
            sh(cmd)
        except KeyboardInterrupt:
            import glob
            for f in glob.glob(s_fasta + ".*"): rm(f)
            raise
    else:
        log.warn("NOT running cmd:\n%s\n because %s.nin is up to date" % (cmd, s_fasta))
    blast_file = ""
    to_query_dir = blast_cfg.get("o", "F").upper() != "F"
    if blast_cfg.get("o", "F").upper() not in ("T", "F"):
        # if it's a file, use the file. otherwise it's a dir, need to
        # create the filename and append it to the dir
        blast_file = blast_cfg["o"]
        if blast_file.endswith(".blast"):
            log.error("using file past in on -o: %s" % blast_cfg["o"])
        else:
            blast_file = get_blast_file(q_fasta, s_fasta, blast_cfg["o"])
            log.error("using directory passed in on -o with exiting file: %s" %
                     blast_file)
    else:
        if full_name:
            blast_file = blast_file.rstrip(".blast") \
                  + "_params__" \
                  + "__".join(["%s_%s" % p for p in sorted(blast_cfg.items())
                               if not p[0] in ("i", "d")]) \
                  + ".blast"
        blast_file = get_blast_file(q_fasta, s_fasta, to_query_dir)

    if blast_log is None:
        blast_log = blast_file + ".log"
    blast_cfg.update({"o": blast_file})
    params = add_dash(blast_cfg)

    params = ["%s %s" % (p, v) for p, v in sorted(params.items())]
    cmd = blastall + " " + " ".join(params)

    if not (is_current_file(blast_file, q_fasta) \
                and is_current_file(blast_file, s_fasta) \
                and is_same_blast_params(blast_file, cmd)):

        fh = open(blast_file + ".cmd", "w")
        fh.write(cmd)
        fh.close()
        try:
            sh(cmd, blast_log=blast_log)
        except:
            rm(blast_file)
            rm(blast_file + ".cmd")
            raise

        if os.path.exists(blast_file):
            lines = sum(1 for line in open(blast_file))
            log.debug("\n\n%s lines of blast output sent to %s" % (lines, blast_file))
        else:
            log.error("\n\nERROR: blast not run")
            if not blast_log is None:
                log.error(open(blast_log).read())
    else:
        log.error("NOT running cmd:\n%s\n because %s is up to date" % (cmd, blast_file))

def check_args(args):
    """
    >>> args = {'i': 'a.fasta', 'd': 'b.fasta'}
    >>> check_args(args)
    >>> args
    {'i': 'a.fasta', 'p': 'blastn', 'd': 'b.fasta', 'a': '4'}
    """
    if not "p" in args: args["p"] = "blastn"
    assert "i" in args, "need to specify a query fasta"
    assert "d" in args, "need to specify a query fasta"
    if not "a" in args: args["a"] = "4"


def handle_temps(args):
    """allow the query (-i) and subject (-d) to be specified as e.g.:
        a.fasta['chr1'][200:500]
    this will create a temporary file of just that chromosome and/or region
    and blast it 
    >>> args = {'i': 'tests/a.fasta', 'd': 'tests/a.fasta'}
    >>> handle_temps(args)
    >>> args['i']
    'tests/a.fasta'

    >>> args['d']
    'tests/a.fasta'

    >>> args['i'] = 'tests/a.fasta["chr2"]'
    >>> try: handle_temps(args) #doctest: +ELLIPSIS
    ... except Exception, e: 
    ...     assert 'no fasta with name' in str(e)

    >>> args['i'] = 'tests/a.fasta[chr1]'
    >>> handle_temps(args) 
    >>> args['i']
    'tests/a.chr1.fasta'
    
    >>> args['d']
    'tests/a.fasta'


    >>> args['d'] = 'tests/a.fasta[chr1][20:25]'
    >>> handle_temps(args) 
    >>> args['d']
    'tests/a.chr1_20_25.fasta'

    >>> open(args['d']).read() #doctest: +NORMALIZE_WHITESPACE
    '>chr1\\nGGGGGG\\n'

    >>> rm(args['d'])
    >>> rm(args['i'])
    """
    def _h(fname):
        if not "[" in fname: return fname
        start = None
        fname = fname.split("[")
        d = os.path.dirname(fname[0])
        if len(fname) == 3:
            fa, seqid, start_stop = [x.rstrip(']').strip("'\"") for x in fname]
            start, stop = [int(x) for x in start_stop.split(":")]
            out_name = os.path.basename(os.path.splitext(fa)[0]) + \
                           (".%s_%s_%s" % (seqid, start, stop)) + ".fasta"
        else:
            fa, seqid = [x.rstrip(']').strip("'\"") for x in fname]
            out_name = os.path.basename(os.path.splitext(fa)[0]) + \
                           (".%s" % (seqid, )) + ".fasta"

        
        assert os.path.exists(fa), fa
        out_name = os.path.join(d, out_name)
        if os.path.exists(out_name) and is_current_file(out_name, fa):
            return out_name

        fh = open(fa, 'rb')
        log.debug('creating sub-file %s' % out_name)
        out = open(out_name, 'wb')
        header = None
        seq = ""
        try:
            for line in fh:
                if header is None:
                    if line[1:].strip() != seqid: continue
                    header = line[1:].strip()
                    print >>out, '>%s' % header
                    continue
                elif line[0] == '>': break
                if start is None:
                    print >>out, line,
                elif len(seq) <= (stop - start):
                    # just hold the entire thing in memory and
                    # snip it at the end.
                    seq += line.rstrip()
            if header is None:
                out.close()
                try:
                    os.unlink(out_name)
                except: pass
                raise Exception("no fasta with name %s containing seq %s"
                                % (fa, seqid))
            if start is not None:
                print >>out, seq[max(0, start - 1):stop]
            fh.close()
            out.close()
            return out_name
        except:
            if os.path.exists(out_name):
                os.unlink(out_name)
            raise

    args['i'] = _h(args['i'])
    args['d'] = _h(args['d'])


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        sh("/usr/bin/blastall")
        print """\
   this script will generally do exactly the same as blastall
   except it will attempt to create the name of the output blast
   file from the input fasta files, and the blast parameters.
   it will also run formatdb with the correct -p parameter based
   on the type of blast requested.
   it also saves a file: a_vs_b.blast.cmd that stores the exact
   blast command used to generate a_vs_b.blast
   in addition, if a.fasta an b.fasta are older than a_vs_b.blast
   and the params have not changed, the blast will not run as
   the current blast file is up to date.

   additional args provided by this script are:
            -o T
         or 
            -o F

         in the former case the blast output file will be created
         from the names of the input fasta files and placed in the
         directory of the query fasta
         in the latter case, the blast file will go to the current
         directory

            --full_name T

         if specified, this will include the blast params in the name
         of the output blast file. e.g.: a__vs_b__params__m_8__W_15.blast
           for -m -8 -W 15
    """
        sys.exit()
    args = dict((sys.argv[i].lstrip("-") , sys.argv[i + 1].rstrip()) \
                          for i in range(1, len(sys.argv), 2))

    try:
        f = args.pop("full_name")
        full_name = not f.lower() in ("f", "0")
    except:
        full_name = False

    if not "i" in args:
        print "need to specify a query fasta (-i)"
        sys.exit()
    if not "d" in args:
        print "need to specify a subject fasta (-d)"
        sys.exit()

    handle_temps(args)
    
    blast(args, full_name=full_name)
