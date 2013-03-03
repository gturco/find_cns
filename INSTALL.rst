CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::

Installation Guide
=================

The easiest way to install the pipeline is to:

  1) Install git and clone this package
  
  + `Install git here <http://git-scm.com/downloads>`_ (click the `link <http://git-scm.com/downloads>`_ and follow instructions for your system )::

       git clone git://github.com/gturco/find_cns.git
       cd find_cns
  
  2) **SEE REQUIREMENTS SECTION**

  3) run bootstrap Installation ::

         ARCHFLAGS='-arch i386 -arch x86_64' python bootstrap.py
  
  + If step 3 fails try installing the failed package (instructions in Manual Installation section) then re-run bootstrap

  4) run test case `read how to run here <https://github.com/gturco/find_cns/blob/master/README.rst#id4>`_::
        
       source ../cns_pipeline/bin/activate
       cd pipeline
       sh run.sh 


Requirements
=============
Install the following (if not already installed) 
**Instructions** are provided below `troubleshooting help available here <http://www.thisisthegreenroom.com/2011/installing-python-numpy-scipy-matplotlib-and-ipython-on-lion/>`_

     - Python 2.7
     - gfortran
     - Blas
     - LAPACK
     - easy_install and pip
     - glpk and scip
     - virtualenv
**Instructions below**

**Special Requirements for Mac**

   + Install homebrew::

      ruby <(curl -fsSkL raw.github.com/mxcl/homebrew/go)

   + Add homebrew to your .bash_profile::

      touch ~/.bash_profile (if this is a new file check ~/.profile or ~/.bashrc)
    Once one of the profiles are open add ``PATH="usr/local/bin/:${PATH}"`` and then ``export PATH`` to a new line
    Then open a new terminal window and type ``brew`` to make sure it was properly installed

  + Install and download `Xcode here <https://itunes.apple.com/us/app/xcode/id497799835?ls=1&mt=12>`_

**Requirements for Mac and Unix**

  + Python version >= 2.7 (you can use `pythonbrew <https://github.com/utahta/pythonbrew/>`_ to install python)::
        
        on mac:
          brew install python --framework --universal
          (edit .bash_profile)
          export PATH=/usr/local/share/python:$PATH
        on ubuntu/debian:
          sudo apt-get update; sudo apt-get install python2.7 python2.7-dev python-setuptools
        on centos/redhat:
          http://toomuchdata.com/2012/06/25/how-to-install-python-2-7-3-on-centos-6-2/

  + gfortran::

        on mac:
            brew install libblas gfortran
        on ubuntu/debian:
            sudo apt-get install build-essential liblas-dev liblapack-dev gfortran
        on centos/redhat:

  + `BLAS <http://www.netlib.org/blas/>`_::
      
        on mac:
          http://pheiter.wordpress.com/2012/09/04/howto-installing-lapack-and-blas-on-mac-os/
        on ubuntu/debian:
          sudo apt-get install libblas-dev
        on centos/redhat:
          sudo yum install blas-devel

  + `LAPACK <http://www.netlib.org/lapack/>_`::

      on mac:
        http://pheiter.wordpress.com/2012/09/04/howto-installing-lapack-and-blas-on-mac-os/
      on ubuntu/debian:
        sudo apt-get install liblapack-dev
      on centos/redhat:
        sudo yum install lapack-devel

  + `GEOS <http://trac.osgeo.org/geos/>`_::

        on mac:
           brew install geos
        on ubuntu/debian:
          sudo apt-get install libgeos-dev
        on centos/redhat:
          sudo yum install geos

    + PIP ::
      
        sudo easy_install pip
 
  + `virtualenv <http://pypi.python.org/pypi/virtualenv/>`_::

        sudo pip install virtualenv
        virtualenv --distribute cns_pipeline --python=python2.7

  + `scip <http://scip.zib.de/download.shtml>`_ Download `here <http://scip.zib.de/download.shtml>`_ choose operating system and **accept user agreement** on next page::
        
        #may need to scp from your computer to server
        unzip scip-x.x.x
        mv scip-x.x.x cns_pipeline/bin/scip

        if on ubuntu/debian need unzip:
          sudo apt-get install unzip (add if not installed)

  + `gpkl <ftp://ftp.gnu.org/gnu/glpk/>`_::
      
        wget glpk-newest_version.tar.gz
        tar -xvzf <somepath>/glpk-newest_version.tar.gz
        cd glpk-newest_version
        ./configure
        make
        sudo make install


      
Manual Installation if bootstrap (step 3) fails
===================================
bootstrap.py runs the  commands below.  If you are having trouble installing one of these packages,  use the links provided.
`troubleshooting numpy and scipy <http://www.thisisthegreenroom.com/2011/installing-python-numpy-scipy-matplotlib-and-ipython-on-lion/>`_


**Python packages**

- First **activate** your virtualenv so everything downloads to your  ``cns_pipeline/bin``::
      
    virtualenv --distribute cns_pipeline --python=python2.7
    (creates folder if not already created)
  
  Then activate::

     source cns_pipeline/bin/activate
     (to deactivate just type: deactivate)

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

- `pandas <http://pandas.pydata.org/>`_::

    pip install pandas

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

-if on mac::

    brew install wget

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

