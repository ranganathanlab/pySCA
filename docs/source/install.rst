============
Installation
============

1. Downloading the code
=======================
The pySCA package, tutorials, and associated scripts are available for download
from the `GitLab repository <https://gitlab.com/sudorook/pySCA>`_.

2. Package Dependencies
=======================
Before running pySCA, you will need to download and install the following
(free) packages:

    1) `Anaconda Scientific Python <https://www.anaconda.com/distribution/>`_
       - this package will install Python, as well as several libraries
       necessary for the operation of pySCA (NumPy, SciPy, IPython, and
       Matplotlib).

    2) Biopython - this can be done in two ways:

       a. From the Anaconda command line, run:

          >>> conda install biopython

       b. Or, download and install from the `Biopython website
          <https://biopython.org/wiki/Download>`_.

    3) Install a robust pairwise alignment program. Either of the below will
       work with pySCA, but in our hands ggsearch is fastest. This is critical
       for the `scaProcessMSA.py` script.

       a. `ggsearch
          <http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml>`_ -
          part of the FASTA software package,

       b. `needle <ftp://emboss.open-bio.org/pub/EMBOSS/>`_ - part of the
          EMBOSS software package,

*The following steps are optional but highly recommended.*	

    4) `PFAM annotations
       <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamseq.gz>`_ -
       the file `pfamseq.txt` contains phylogenetic annotations for PFAM
       sequences. This is necessary if you would like to annotate PFAM
       alignments with taxonomic/phylogenetic information using the
       annotateMSA.py script provided by pySCA. This file is quite large (~10
       GB) and is not included here, but it is available from the PFAM FTP
       site in compressed (\*.gz) format.

    5) `PyMol <https://pymol.org/2/>`_ - necessary if you would like to
       use pySCA's automated structure mapping scripts, and useful for mapping
       the sectors to structure in general.

    6) `mpld3 <http://mpld3.github.io/>`_ - a package that allows more
       interactive plot visualization in Jupyter notebooks. If you choose not to
       install this (optional) package, you will need to comment out the
       `import mpld3` lines at the beginning of the tutorials.


3. Path Modifications
=====================

Following the successful installation of these packages, edit the following
variables in the "PATHS" of `settings.py` to reflect the locations of these
files on your computer. (Change the values of each variable to match however
the files ar found on your system if they differ from the defaults.)

:path2pfamseq: location of the pfam.seq database file (default: `data/pfamseq.txt`)

:path2structures: location of your PDB structures for analysis (default: `data/`)

:path2pymol: location of your PyMOL executable (default: `/usr/bin/pymol`).  If
             you don't know where it is, try running (if on a Unix system):

>>> whereis pymol # requires Unix shell

:path2needle: location of the needle executable, if you have installed EMBOSS
              needle (default: `/usr/bin/`)

:path2output: location of output directory (default: `output/`)


4. Getting started and Running the tutorials
============================================
The `"getting started"`_ section of this documentation provides instructions on
how to run some initial calculations and the tutorials. The basic idea behind
the pySCA code is that the core calculations are performed using a series of
executable Python scripts, and then the results can be loaded and
analyzed/visualized using an Jupyter notebook (or alternatively, Matlab).

All of the tutorials are written provided as Jupyter notebooks. For more on
how Jupyter notebooks work, see: http://ipython.org/notebook.html. Prior to
running the notebook tutorials, you'll need to run the core calculation scripts
that generate the input for the notebooks. One way to do this is with the shell
script "runAllNBCalcs.sh", and there is more information on this in the
`"getting started"`_ section. Once the calculations are completed, you can
begin the tutorial in interactive Python from the command line, by typing:

>>> ipython notebook <NOTEBOOK_NAME_HERE>

.. _"getting started": get_started.html
