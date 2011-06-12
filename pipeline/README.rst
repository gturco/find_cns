CNS Pipeline
============

:Author: Brent Pedersen (brentp)
:Email: bpederse@gmail.com
:License: MIT

.. contents ::

Requirements
============

 + `flatfeature <http://github.com/brentp/flatfeature/>`_
   (check with git and run ``sudo python setup.py install``)

 + `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (download latest .tar.gz; configure; make; make install) and adjust path in quota.sh)

 + `quota-align <http://github.com/tanghaibao/quota-alignment>`_
   (checkout with git and adjust path in quota.sh)

 + `numpy <http://github.com/numpy/numpy/>`_ checkout and run ``sudo python setup.py install``

 + `pyfasta` (``sudo easy_install -UZ pyfasta`` you will have latest from pypi).

 + `shapely` (``sudo apt-get install libgeos-dev``, ``sudo easy_install -UZ 'shapely==1.0.0'``)

 + `processing` (``sudo easy_install -UZ processing``)

 + `bl2seq`

 + `blastall`

 + `bpbio`



Run
===

 + *Once only*: edit quota.sh to correct path for ``quota-alignment``
 + for each organism, use export_to_bed.pl to get data out of CoGe. e.g.::

    perl scripts/export_to_bed.pl \
                          -fasta_name rice_v6.fasta \
                          -dsg 8163 \
                          -name_re "^Os\d\dg\d{5}$" > rice_v6.bed

   where ``dsg`` is from CoGe OrganismView and the prefix for the .bed and
   .fasta file **must be the same** (in this case ``rice_v6``).
   You likely need to run this on new synteny and then copy the .bed and
   .fasta files to the ``data/`` directory.
   The -name_re regular expression is not required, but in this case, it will
   prefer the readable Os01g101010 names over the names like m103430.

 + edit quota.sh to the correct `ORGA`, `ORGB`, `QUOTA`
 + run `sh run.sh` # that will call quota.sh (this will take a long time as it's doing
   a full blast (lastz) and then all of quota align, then cns pipeline).
 + this will create png's for the dotplots. check those to make sure the quota-blocks look correct.
