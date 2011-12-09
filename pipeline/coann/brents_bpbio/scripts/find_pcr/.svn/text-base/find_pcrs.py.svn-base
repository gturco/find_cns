"""
given a fasta file with alternating query, subject cns-sequences. e.g:
    >q__Bd1|59025115|59025142|1|55266333|55266306
    AAAGGAAAACCAACTTCAGCAATCAGAA
    >s__Bd1|59025115|59025142|1|55266333|55266306
    TTCTGATGGCTGAAGTTGGTTTTCTTTT

report pairs which will make good pcr sites.

"""
import sys
import os
import os.path as op
import subprocess
import string
sys.path.insert(0, op.dirname(__file__))


commandl = \
""" echo 'PRIMER_TASK=pick_left_only\nSEQUENCE=%s\n=' | /opt/src/primer3/primer3-1.1.3/src/primer3_core | egrep 'PRIMER_LEFT\\b|PRIMER_LEFT_SEQUENCE'"""
commandr = \
""" echo 'PRIMER_TASK=pick_right_only\nSEQUENCE=%s\n=' | /opt/src/primer3/primer3-1.1.3/src/primer3_core | egrep 'PRIMER_RIGHT\\b|PRIMER_RIGHT_SEQUENCE'"""

fasta = sys.argv[1]
outdir = op.dirnname(fasta)

cns_fh = open(fasta, "r")
out_file = sys.stdout
line = True
goods = []
bads = []

_trans = string.maketrans("ATCGNatcg", "TAGCNtagc")
rcomp = lambda seq: seq.translate(_trans)[::-1]

seen = {}
i = 0
while 1:
    headerq = cns_fh.readline()[1:].strip()
    seqq    = cns_fh.readline().strip().upper()
    headers = cns_fh.readline()[1:].strip()
    seqs    = cns_fh.readline().strip().upper()
    if not seqs: break
    if not (headerq.startswith("q__") and headers.startswith("s__")):
        print headerq, headers
        raise

    assert set(seqs.lower()).issubset("atcgnx"), (seqs.lower(), headers)
    assert set(seqq.lower()).issubset("atcgnx"), (seqq.lower(), headerq)
    #print >>sys.stderr, headerq, headers, seqq, seqs
    if (headerq, headers) in seen: continue
    seen[(headerq, headers)] = 1

    #print command % seq
    # save all the primers for both right and left.
    primer_list = []
    # to allow sorting based on A/B/C hits.
    ABCs = []
    for lr, command in (('left', commandl), ('right', commandr)):
        pq = subprocess.Popen(command % seqq, shell=True, stdout=subprocess.PIPE)

        ps = subprocess.Popen(command % seqs, shell=True, stdout=subprocess.PIPE)
        primerq = pq.communicate()[0].rstrip()
        primers = ps.communicate()[0].rstrip()

        # neither have a result.
        if not (primers or primerq): continue

        # just get the sequence.
        tailq, tails = None, None
        qpos, spos = None, None
        if primerq:
            a = primerq.split()
            #print >>sys.stderr, a
            qpos = map(int, a[1][a[1].find('=') + 1:].split(","))
            primerq =  a[0][a[0].find('=') + 1:]
            if lr == 'left':
                tailq = primerq[-5:]
            else:
                tailq = rcomp(primerq[:5])

        if primers:
            a = primers.split()
            #print >>sys.stderr, a
            spos = map(int, a[1][a[1].find('=') + 1:].split(","))
            primers =  a[0][a[0].find('=') + 1:]
            if lr == 'left':
                tails = primers[-5:]
            else:
                tails = rcomp(primers[:5])

        
        quality = 'C'   
        if primerq == primers:
            quality = 'A'
            primer_list.append(quality + ':' + lr + ":" + primerq)
            ABCs.append(quality)
        else:
            if tailq:
                posq = seqq.find(tailq)
                assert posq != -1, (tails, seqs, command % seqs)
                poss = seqs.find(tailq)
                if poss != -1 and abs(posq - poss) < 2:
                    quality = 'B'
                primer_list.append(quality + ':' + lr + ":" + primerq)
                ABCs.append(quality)
            if tails:
                poss = seqs.find(tails)
                assert poss != -1, (tails, seqs, command % seqs)
                posq = seqq.find(tails)
                if posq != -1 and abs(posq - poss) < 2:
                    quality = 'B'
                primer_list.append(quality + ':' + lr + ":" + primers)
                ABCs.append(quality)
        

    if primer_list != []:
        key = headerq[3:]
        #current = shelf[key]
        #abcs = "|".join(sorted(set(ABCs)))
        #current['pcr'] = {'primer_list': primer_list, 'quality': abcs }
        #shelf[key] = current
        assert key == headers[3:]

        print >>out_file, key, "|".join(sorted(set(ABCs))), seqq, seqs, "|".join(primer_list)

