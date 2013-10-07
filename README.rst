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

Cite
============

Turco, G., Schnable, J. C., Pedersen, B., & Freeling, M. Automated conserved noncoding sequence (CNS) discovery reveals differences in gene content and promoter evolution among grasses. Frontiers in Plant Science, 4, 170. `Link to Plant CNS Paper <http://www.frontiersin.org/plant_genetics_and_genomics/10.3389/fpls.2013.00170/abstract>`_


Datasets
============
CNS Datasets Available `here <http://figshare.com/articles/CNS_Discovery_Pipeline_Supporting_Data/107054>`_ for:

  - Automated thaliana_v10 thaliana_v10 cns
  - Automated rice maize cns 
  - Automated rice sorghum cns
  - Automated sorghum sorghum
  - Automated rice setaria cns
  - Automated setaria setaria

Please use citation above 

Installation
============

  - Everything is predownloaded except scip in the `iplant atmosphere <https://atmo.iplantcollaborative.org/login/>`_ image find_cns_pipeline (emi-8EF728EB)
  - Download the most recent code here::
      
      git clone https://github.com/gturco/find_cns.git

  - run bootstrap code and add scip to cns_pipline/bin/::

       python bootstrap.py

**Required Dependencies** 

  - Read `INSTALL file <https://github.com/gturco/find_cns/blob/master/INSTALL.rst>`_ for instructions

Running the Pipeline
====================

**Obtaining Input Files**

 - Download Fasta and gff from `CoGe OrganismView <http://genomevolution.org/CoGe/OrganismView.pl>`_ for each organism 
 - Make sure to **UNCLICK**  "Do not generate features for ncRNA genes (CDS genes only)" when downloading gff

 Convert gff to `Bed format <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_::

      python scripts/gff_to_bed.py --re "^Os\d\dg\d{5}" --gff rice_v6.gff  --fasta 9109.faa --out rice_v6

 - The -re regular expression is not required, but in this case, it will prefer the readable Os01g101010 names over the names like m103430.
 - the --out is the root word for the Fasta and Bed outfiles since they **must be the same name** (in this case ``rice_v6.fasta and rice_v6.bed``)
 - Fasta File (it is recommended to run `50x mask repeat <http://code.google.com/p/bpbio/source/browse/trunk/scripts/mask_genome/mask_genome.py>`_)


**Runing Pipeline**


 - Create a new directory (under data) name for organisms being compared ORGA_ORGB  eg ``mkdir data/rice_v6_sorghum_v1``
 - Add Fasta and Bed files to directory
 - Edit run.sh file ``touch run.sh`` change `ORGA`, `ORGB`, `QUOTA`, `SDGID`, `QDGID` #DGID found in CoGe
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

