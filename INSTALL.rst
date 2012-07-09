CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::

Installation Guide
=================

The easiest way to install the pipeline is to:
  
  1) git clone this package ::

       git clone : git://github.com/gturco/find_cns.git

  2) Download the requirements below
  3) run Installation ::

         ARCHFLAGS='-arch i386 -arch x86_64' python bootstrap.py
  
  + If step 3 fails then try installing the failed package using one of the steps below and re-doing step 3


Requirements
===========
The fllowing are absolutely requires otherwise installation WILL fail

  + Python version >= 2.7 (you can use `pythonbrew <https://github.com/utahta/pythonbrew/>`_ to install python)
  + `git <https://help.github.com/articles/set-up-git>`_
  + gfortran (type gfortran to check if on machine or install with `homebrew <https://github.com/mxcl/homebrew/wiki/Installation>`_ or `xcode <https://developer.apple.com/xcode/>`_)
  + `GEOS <http://trac.osgeo.org/geos/>`_

  IF not using bootstrap.py:

  + `activate virtualenv <http://pypi.python.org/pypi/virtualenv/>`_::

      python virtualenv.py --distribute cns_pipeline --python=Python2.7
      source cns_pipeline/bin/activate


Installation
============
To use the CNS pipeline the following python packages and c packages are required

**Python packages**

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


**C packages**

- `(NON-blast+) blast <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/>`_
   download latest blast from  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/::

    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.5/blast-2.2.5-ia32-linux.tar.gz
    tar -xvzf <somepath>/blast-X.X.X-XXXX.tar.gz
    mv <somepath>/blast-XX.X.X/ cns_pipeline/bin/ #(change path in run.sh file if diff)

- `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (`install instructions <http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#install>`_ ) and adjust path in quota.sh)::

    wget http://www.bx.psu.edu/~rsharris/lastz/newer/lastz-1.03.02.tar.gz
    tar -xvzf <somepath>/lastz-distribute-X.XX.XX.tar.gz
    cd <somepath>/lastz-distrib-X.XX.XX/src
    make
    LASTZ_INSTALL=/usr/local/bin/ make install

