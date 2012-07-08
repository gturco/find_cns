CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::

Requirements
===========
  + Python version >= 2.7
  + PIP
  + libgeos_c
  + gfortran
  + git
  + `GEOS <http://trac.osgeo.org/geos/>`_

.. image:: http://upload.wikimedia.org/wikipedia/commons/0/08/Pipeline_git.png

Python Installation
============

- `numpy <http://www.scipy.org/Download/>`_::

    pip install numpy

- `processing <http://pypi.python.org/pypi/processing/>`_::

    pip install processing

- `shapely <http://toblerity.github.com/shapely/manual.html>`_::

    pip install shapely

- `pyfasta <http://pypi.python.org/pypi/pyfasta/>`_::

    pip install pyfasta

- `scipy <http://www.scipy.org/Installing_SciPy/>`_::

    pip install scipy

- `Cython <http://www.cython.org/#download>`_::

    pip install Cython

- `flatfeature <https://github.com/brentp/flatfeature.git>`_::

    pip install git+https://github.com/brentp/flatfeature.git

- `quota-align <https://github.com/tanghaibao/quota-alignment>`_::
  
    git clone https://github.com/tanghaibao/quota-alignment.git 
    mv quota-alignment  cns_pipeline/bin/
    (change path in quota.sh if not moved to cns_pipeline/bin)



- `gffparser <https://github.com/chapmanb/bcbb/tree/master/gff>`_::

    git clone https://github.com/chapmanb/bcbb.git
    cd gff
    python setup.py install

- `bpbio <http://code.google.com/p/bpbio/>`_::

    cd pipeline/coann/brents_bpbio/biostuff/
    python setup.py install
    cd pipeline/coann/brents_bpbio/blasttools/blast_misc/
    python setup.py install
    cd pipeline/coann/brents_bpbio/biostuff/co-anno/
    python setup.py install

- test


C Installation
============

 + `blast <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>`_
   (download latest at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/  and run)

 + `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (download latest .tar.gz; configure; make; make install) and adjust path in quota.sh)
