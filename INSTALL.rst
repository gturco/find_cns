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

.. image:: http://upload.wikimedia.org/wikipedia/commons/0/08/Pipeline_git.png

Python Installation
============
  + numpy : ``pip install numpy``
  + processing: ``pip install processing``
  + shapely: ``pip install shapely``
  + pyfasta: ``pip install pyfasta``
  + scipy: ``pip install scipy``
  + Cython: ``pip install Cython``
  + flatfeature: ``pip install git+https://github.com/brentp/flatfeature.git``
  + quota-align: ``git clone https://github.com/tanghaibao/quota-alignment.git`` (mv to /cns_pipeline/bin/)  (checkout with git and adjust path in quota.sh)
  + gffparser: ``git clone https://github.com/chapmanb/bcbb.git``
    ``cd gff``
    ``python setup.py install``
  + ``cd pipeline/coann/brents_bpbio/biostuff/ \n python setup.py install``
  + ``cd pipeline/coann/brents_bpbio/blasttools/blast_misc/ \n python setup.py install``
  + ``cd pipeline/coann/brents_bpbio/biostuff/co-anno/ \n python setup.py install``




C Installation
============

 + `blast <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>`_
   (download latest at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/  and run)

 + `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (download latest .tar.gz; configure; make; make install) and adjust path in quota.sh)
