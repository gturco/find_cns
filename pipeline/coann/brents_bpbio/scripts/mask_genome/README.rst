==============
Genome Masking
==============

given a full genome self-`blast` and a `fasta` file, create
a new masked genome where any basepair occuring more than
`cutoff` times is masked. the new file will be written
in the same directory as the original fasta with a name like
"rice.masked.50.fasta" where rice.fasta was the input file
and 50 the cutoff.

Blast
=====
typical parameters used for the blast are
::

  $ bblast.py -i $FASTA -d $FASTA -a 4 -e 0.001 -m 8 -o ${ORG}_${ORG}.blast

Soft Masking
============
by default the genome is 'hard' masked with 'X', that is any repetitive basepair
is masked with 'X'. it is possible to specify a soft mask with "-m SOFT " on the
commandline. this will convert each sequence to upper case and then mask repetetive
values with the lower-case equivalent. tools including BLAST can be told to ignore the lower-cased basepairs ("-U T" on the command-line for blastall).

Example
=======
an entire mask. starting from just the genomic fasta can be run with
::

    ORG=thaliana_v7
    FA=${ORG}.fasta
    BL=${ORG}_${ORG}.blast
    bblast.py -p blastn -i $FA -d $FA -m 8 -a 8 -o $BL
    python mask_genome.py -b $BL -f $FA -o $ORG -c 50


after running the self-self blast on the thaliana_v7.fasta, this will create a
new fasta file "thaliana_v7.masked.50.fasta" with all basepairs covered by 
more than 50 blast hits masked to 'X'. 
An hdf5 file named 'copy_count.h5' will be created (or added to) with a group
of 'thaliana_v7' which contains an array for each of the 5 chromosomes. Each array
contains as many entries as there are basepairs in the corresponding chromosome
and each entry corresponds to the number of blast hits covering that basepair.

Plotting
========
the script mask_plot.wsgi is included. it assumes that modwsgi
is installed for the apache webserver (but can also be adapted to
cgi or mod-python).
given an apache conf like
::
    
    $ WSGIScriptAlias /copy_count /full/path/to/mask_genome/mask_plot.wsgi

and `thaliana_v7` run through the masking and the resulting copy_count.h5 
placed in the same directory as mask_plot.wsgi, a url request like
::

    http://localhost/copy_count/?org=thaliana_v7&xmin=12&xmax=90000&width=800&seqid=2

will generate an image 800 pixels wide of chromosome `2` from basepair `12`
to `90000` the y-range of the plot is from 0 to 200 (easily changed) and
the line drawn is the number of blast-hits along the chromosome at each
basepair location. a horizontal red line is drawn at 50 basepairs as that is
a commonly used value for masking.

Requirements
============

 * hd5 c library and headers
 * pytables (easy_install - able)
 * pyfasta (easy_install - able)
 * numexpr (easy_install - able)

