============
Installation
============

1. Downloading the code
=======================

The pySCA package, tutorials, and associated scripts are available for download
from the `GitHub repository <https://github.com/ranganathanlab/pySCA>`_.


2. Package Dependencies
=======================

Before running pySCA, you will need to download and install the following
(free) packages:

    1) `Anaconda Scientific Python <https://www.anaconda.com/distribution/>`_
       - this package will install Python, as well as several libraries
       necessary for the operation of pySCA (NumPy, SciPy, Jupyter, and
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
       <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamseq.txt.gz>`_ -
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

Your requirements will vary depending on the size of your sequence alignments,
but as a rule of thumb, the toolbox is best used on a system with at least 8 GB
of RAM. Less will run the risk of processes using more memory than available
and subsequently getting killed by the operating system's scheduler.


3. Path Modifications
=====================

Following the successful installation of these packages, you may edit some
variables for your convenience before installing the pySCA package. The
following variables in the "PATHS" of `settings.py` to reflect the locations of
these files on your computer. (Change the values of each variable to match
however the files ar found on your system if they differ from the defaults.)

:path2pfamseq: location of the pfamseq text file (default: `data/pfamseq.txt`)

:path2pfamdb: location of the pfamseq database (default: `data/pfamseq.db`)

:path2structures: location of your PDB structures for analysis (default:
                  `data/`)

:path2pymol: location of your PyMOL executable (default: `/usr/bin/pymol`). If
             you don't know where it is, try running (if on a Unix system):

.. code-block:: shell
 
  $ whereis pymol

:path2needle: location of the needle executable, if you have installed EMBOSS
              needle (default: `/usr/bin/`)

:path2output: location of output directory (default: `output/`)

The functions that use these global settings can also handle them as
command-line arguments. This file is present for convenience.

Note: the :code:`$` at the start of each line typically denotes that the
supplied code is to be executed in terminal, and similarly :code:`>>>` denotes
code for the Python REPL. Don't enter the :code:`$`  or :code:`>>>` when
copying and running each line.


4. Installation
===============

The analysis scripts found in the bin/ directory and the SCA toolbox in pysca/
can now be installed. To install them system-wide, from the base of the
repository::

  $ sudo pip install .

Or, if `pip` is not installed, run::

  $ sudo python setup.py install

If you'd rather install this locally (i.e. without superuser permissions), you
can instad run::
  
  $ pip install --user .

Or::

  $ python setup.py install --user

Note that to use locally installed scripts, the installation directory needs to
be in the system PATH. To check whether that is the case, run::

  $ echo $PATH | grep --color=auto "$(python -m site --user-base)/bin"

If the installation directory is highlighted in the output, then the PATH is
configured correctly. If it is not found, then it needs to be added manually.
Open you shell configuration file (e.g. .bashrc) and add the directory to the
PATH variable by appending the following line::

  export PATH="$PATH:$HOME/.local/bin"

The exact path (the text following the semicolon) may differ on your system,
but it can easily be found by running `echo $(python -m site --user-base)/bin`.

Now, with the pySCA code installed, each of the commands found in bin/ can now
be run from the command line.


4. Getting Started and Running the Tutorials
============================================

The :doc:`"getting started" <get_started>` section of this documentation
provides instructions on how to run some initial calculations and the
tutorials. The basic idea behind the pySCA code is that the core calculations
are performed using a series of executable Python scripts, and then the results
can be loaded and analyzed/visualized using an Jupyter notebook (or
alternatively, MATLAB).

All of the tutorials are written provided as Jupyter notebooks. For more on
how Jupyter notebooks work, see: `<https://jupyter.org>`_. Prior to running the
notebook tutorials, you'll need to run the core calculation scripts that
generate the input for the notebooks. One way to do this is with the shell
script "runAllNBCalcs.sh", and there is more information on this in the
:doc:`"getting started" <get_started>` section. Once the calculations are
completed, you can begin the tutorial in interactive Python from the command
line, by typing:

>>> jupyter notebook <NOTEBOOK_NAME_HERE>
