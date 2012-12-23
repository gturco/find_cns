from biostuff import BlastLine, BlastFile

some_attrs = ('qstart', 'qstop', 'sstart', 'sstop', 'pctid', 'score', 'query',
        'subject')

def test_blastfile():
    f = "tests/data/tabd.blast" 
    bf = BlastFile(f)
    fh = open(f, 'r')

    # iterate via python and c and check each line is the same.
    for line, b in zip(fh, bf):
        bl = BlastLine(line)
        assert isinstance(b, BlastLine)
        assert bl == b

    i = 0
    for c in bf:
        i += 1
        assert isinstance(c, BlastLine)
    assert i == len(open(f).readlines())

    del bf

def test_blastfile_list():
    f = "tests/data/tabd.blast" 
    blasts = list(BlastFile(f))
    assert len(blasts) == len(open(f).readlines())

def test_blastline():
    f = "tests/data/tabd.blast" 
    blasts = []
    for line in open(f):
        b = BlastLine(line)
        blasts.append(BlastLine(line))

    yield check_type, blasts, ('qstart', 'qstop', 'sstart', 'sstop',
        'nmismatch', 'ngaps'), int

    yield check_type, blasts, ('evalue', 'score', 'pctid'), float
    yield check_type, blasts, ('query', 'subject'), str
    

def check_type(blasts, attrs, klass):
    for b in blasts:
        for attr in attrs:
            assert isinstance(getattr(b, attr), klass)


def test_query_subject_props():
    f = "tests/data/tabd.blast" 
    line = BlastLine(open(f).readline())
    line.query = "asdf"
    line.subject = "dddd"
    assert line.query == "asdf"
    assert line.subject == "dddd"
    assert "asdf" in line.to_blast_line()
    assert "dddd" in line.to_blast_line()

def test_to_string():
    f = "tests/data/tabd.blast" 
    for line in open(f):
        a = BlastLine(line)
        b = BlastLine(a.to_blast_line())

        # works better than string comparison because of floats.
        for attr in some_attrs:
            assert getattr(a, attr) == getattr(b, attr), (a, b, attr)

def test_pickle():
    import cPickle
    f = "tests/data/tabd.blast" 
    line = BlastLine(open(f).readline())
    
    d = cPickle.dumps(line, -1)

    loaded = cPickle.loads(d)

    for k in BlastLine.attrs:
        assert getattr(loaded, k) == getattr(line, k)
    loaded.query = "asdf"

    assert loaded.query != line.query
