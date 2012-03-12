CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::


Description
===========
Python application to automate large-scale identification of conserved noncoding sequences (CNS) between two `usefully diverged <http://genomevolution.org/wiki/index.php/Useful_divergence>`_ plant species.
This application works by first attempting to correct annotation errors between the two species. It then condenses local duplicates and finds syntenic regions based on `ploidy relationships <https://github.com/tanghaibao/quota-alignment>`_. BLAST is then applied to the syntenic regions between the two species to find CNSs. CNSs are found through blastn at an e-value less than or equal to a 15/15 exact base pair match. Nonsyntenic CNSs are removed along with CNS with hits to known RNA or exons.

.. image:: http://genomevolution.org/wiki/images/6/6e/Peach-Chocolate-example.png

Installation
============
 + Python version >= 2.7

 + `blast <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>`_
   (download latest and run)

 + `lastz <http://www.bx.psu.edu/~rsharris/lastz/newer/>`_
   (download latest .tar.gz; configure; make; make install) and adjust path in quota.sh)

 + `quota-align <http://github.com/tanghaibao/quota-alignment>`_
   (checkout with git and adjust path in quota.sh)

 + `flatfeature <http://github.com/brentp/flatfeature/>`_
   (check with git and run ``sudo python setup.py install``)

 + `pyfasta <https://github.com/brentp/pyfasta>`_ (``sudo easy_install -UZ pyfasta`` you will have latest from pypi).

 + `shapely <http://pypi.python.org/pypi/Shapely#downloads>`_ (``sudo apt-get install libgeos-dev``, ``sudo easy_install -UZ 'shapely==1.0.0'``)

 + `numpy <http://github.com/numpy/numpy/>`_ checkout and run ``sudo python setup.py install``

 + `processing <http://pypi.python.org/pypi/processing#downloads>`_ (``sudo easy_install -UZ processing``)

 + `bpbio <http://pypi.python.org/pypi/processing#downloads>`_ (``svn checkout http://bpbio.googlecode.com/svn/trunk/ bpbio-read-only``)
   (run biostuff,coanno and bblast ``sudo python setup.py install``)


Run
===
Inputs
-------

 + Fasta File (it is recommended to run `50x mask repeat <http://code.google.com/p/bpbio/source/browse/trunk/scripts/mask_genome/mask_genome.py>`_)
 + Bed File (supports `UCSC bed format <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_)
 + Converting GFF to Bed
  ``BCBio`` module required::
      
      python scripts/gff_to_bed.py rice_v6.gff >rice_v6.bed


 + If you have access to Coge the fasta and bed file for each organism can be obtained using export_to_bed.pl e.g.::

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


Editing Run File
::::::::::::::::

 + *Once only*: edit quota.sh to correct path for ``quota-alignment``
 + edit quota.sh to the correct `ORGA`, `ORGB`, `QUOTA`
 + under the data directory create a new directory titled ORGA_ORGB with their corresponding beds and fasta eg ``mkdir data/rice_v6_sorghum_v1``
 + run cmd: ``sh run.sh`` #that will call quota.sh (this will take a long time as it's doing a full blast (lastz) and then all of quota align, then cns pipeline).
 + this will create png's for the dotplots. check those to make sure the quota-blocks look correct.

Output files
::::::::::::

 + Query and subject CNS position
 + Missing Exons from ORGA ORGB blast
 + CNS blast to  RNA file
 + CNS blast to proteins file
 + CNS assigned to nearest Ortholog
