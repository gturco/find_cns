CNS Pipeline
============
:Author: Gina Turco (`gturco <https://github.com/gturco>`_), Brent Pedersen (`brentp <http://github.com/brentp>`_)
:Email: gturco88@gmail.com
:License: MIT
.. contents ::

Description
===========
Python application to automate large-scale identification of conserved noncoding sequences (CNS) between two `usefully diverged <http://genomevolution.org/wiki/index.php/Useful_divergence>`_ plant species.
This application works by first attempting to correct annotation errors between the two species using `co-anno <https://github.com/gturco/co-anno>`_. It then condenses local duplicates and finds syntenic regions based on ploidy relationships using `quota-alignment <https://github.com/tanghaibao/quota-alignment>`_. BLAST is then applied to the `syntenic regions <http://genomevolution.org/wiki/index.php/Syntenic_regions>`_ between the two species to find CNSs. CNSs are found through blastn at an e-value less than or equal to a 15/15 exact base pair match (`Kaplinsky et al. <http://www.pnas.org/content/99/9/6147.long>`_). Nonsyntenic CNSs are removed along with CNS with hits to known RNA or exons.
Created in the `Freeling Lab <http://microscopy.berkeley.edu/~freeling/>`_ at UC Berkeley

.. image:: http://upload.wikimedia.org/wikipedia/commons/0/08/Pipeline_git.png

Installation
============

  - Download most recent code here::
      
      git clone https://github.com/gturco/find_cns.git

**Required Dependencies** 

  - Read `INSTALL file <https://github.com/gturco/find_cns/blob/master/INSTALL.rst>`_ for instructions

Running the Pipeline
====================

**Obtaining Input Files**

 - Requires a fasta file and `Bed file <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_ for each organism being compared 
 - Download fasta and gff from `CoGe OrganismView <http://genomevolution.org/CoGe/OrganismView.pl>`_ for each organism 
 - Make sure to **UNCLICK**  "Do not generate features for ncRNA genes (CDS genes only)" 
 - Convert files to correct format::
      python scripts/gff_to_bed.py -re "^Os\d\dg\d{5}" --gff rice_v6.gff --bed rice_v6.bed
      perl -pi -e "s/>gi\|/>/gi" 9109.faa
  
 - The -re regular expression is not required, but in this case, it will prefer the readable Os01g101010 names over the names like m103430.
 - Fasta File (it is recommended to run `50x mask repeat <http://code.google.com/p/bpbio/source/browse/trunk/scripts/mask_genome/mask_genome.py>`_)
 - **RENAME** Fasta file and Bed file **must be the same name** (in this case ``rice_v6``).
 - Move fasta and Bed file to the ``data/`` directory


**Runing Pipeline**


 - edit run.sh to the correct `ORGA`, `ORGB`, `QUOTA`, `SDGID`, `QDSGID`
 - under the data directory create a new directory titled ORGA_ORGB with their corresponding beds and fasta eg ``mkdir data/rice_v6_sorghum_v1``
 - activate screen ``screen``
 - activate virtualenv: ``source ../cns_pipeline/bin/activate``
 - run cmd: ``sh run.sh`` #that will call quota.sh (this will take a long time as it's doing a full blast (lastz) and then all of quota align, then cns pipeline).
 - this will create png's for the dotplots. check those to make sure the quota-blocks look correct.
 - when finshed deactivate virtualenv ``deactivate``

**Output Files**


 - CNSlist (contains start,stop,chr,sequences and 5 prime 3 prime information for gene)
 - Genelist  (one for each ORG, contains the gene start,stop,local dupinfo,orthos,number of CNSs)
 - new genes are named ORG_chr_start_stop in the Genelist
 - CNS with hits to rna or protein are also renamed in the Genelist

