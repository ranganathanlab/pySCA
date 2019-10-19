============
Installation
============

1. Dependencies
===============

Before installing pySCA, you will need to download and install the following
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

**Note**: The :code:`$` at the start of each line typically denotes that the
supplied code is to be executed in terminal, and similarly :code:`>>>` denotes
code for the Python REPL. Don't enter the :code:`$`  or :code:`>>>` when
copying and running each line.

**Important:** The :code:`ggearch36`, :code:`needle`, and :code:`pymol`
programs need to be on the system PATH.

To view your system PATH, run::

  $ echo $PATH

To add directories containing the required prorams to your system path, you
will need to edit your shell configuration file (e.g. `.bashrc` or
`.bash_profile`) found at the base of your user directory. To add a directory
to the system PATH, open up the file and apped the line::

  export PATH="$PATH:<path to directory>"

where `<path to directory>` is replaced with the path to the directory
containing a program you wish to add (e.g. `~/.local/bin`). After saving the
changes, new terminals will use the updated PATH.

**Important:** To add an already-installed program is to the PATH, run::

  $ whereis <program>

to find where `<program>` (e.g. :code:`pymol`) is located, and add its
directory to the system PATH in the manner described above. 

**Important:** Your requirements will vary depending on the size of your
sequence alignments, but as a rule of thumb, the toolbox is best used on a
system with at least 8 GB of RAM. pySCA may run with Less, but there will be a
greater risk when using modestly-sized multeiple sequence alignments of
processes using more memory than available and subsequently getting killed by
the operating system's scheduler.


2. Downloading the Code
=======================

The pySCA package, tutorials, and associated scripts are available for download
from the `GitHub repository <https://github.com/ranganathanlab/pySCA>`_. There
are several options for doing so.

A. Use Git
-----------

If you have :code:`git` installed on your system, you can use it to clone the
repository from GitHub. Run::

  $ git clone https://github.com/ranganathanlab/pySCA.git

For development and troubleshooting purposes, using Git is preferred.

B. Download from the website
----------------------------

If you don't use :code:`git`, you can download the source code from the GitHub
website. Click the green "Clone or download" tab pictured below to obtain the
latest code.

.. image:: _static/github-download-screenshot.png

In the event that you need older versions of the code, you can use the
`releases <https://github.com/ranganathanlab/pySCA/releases>`_ tab on the
GitHub page to download older tagged version.


3. (OPTIONAL) Path Modifications
================================

Before installing pySCA, for your convenience, you may specify default paths in
te `settings.py` file for you to use instead of needed to type them out in the
command line. This part is optional, as every path that you can set in this
file can be specified in the command line.

The following variables in the "PATHS" of `settings.py` to reflect the
locations of these files on your computer.

:path2pfamseq: location of the pfamseq text file (default: `data/pfamseq.txt`)

:path2pfamdb: location of the pfamseq database (default: `data/pfamseq.db`)

:path2structures: location of your PDB structures for analysis (default:
                  `data/`)

:path2output: location of output directory (default: `output/`)


4. Installation
===============

The analysis scripts found in the bin/ directory and the SCA toolbox in pysca/
can now be installed. To install them system-wide, from the base of the
repository::

  $ pip install .

You will of course need to have :code:`pip` installed.

In the event you run in to permissions errors, two options are to either:

A. Install pySCA Locally
------------------------

To install pySCA in your user directory (and without root priviledges), run::

  $ pip install --user .

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

B. Install pySCA System-wide
----------------------------

To install pySCA system-wide, run (as root/administrator)::

   $ sudo pip install .

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
