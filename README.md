# pySCA

![Website Build Status](https://gitlab.com/ranganathanlab/pySCA/badges/master/pipeline.svg)

> Version 7.0
>
> Copyright (C) 2019 Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds
>
> This program is free software distributed under the BSD 3-clause license,
> please see the file LICENSE for details.

pySCA is a Python 3 implementation of the Statistical Coupling Analysis (SCA) for studying the architecture of evolutionary conservation in proteins. The analysis starts with a multiple sequence alignment of a protein family and produces a description of both positional conservation and collectively evolving groups of amino acids (sectors).

## Quick Start

1. **Install pySCA:**
   ```bash
   pip install -e .
   ```

2. **Process an alignment:**
   ```bash
   sca-process-msa alignment.fasta -s 1XYZ --chainID A
   ```

3. **Run SCA calculations:**
   ```bash
   sca-core Outputs/alignment.db.gz --do-sector-id
   ```

For detailed installation instructions, see [INSTALLATION.md](INSTALLATION.md).

For usage instructions, see [USAGE_INSTRUCTIONS.md](USAGE_INSTRUCTIONS.md).

For more information, visit:
- [pySCA Website](https://ranganathanlab.gitlab.io/pySCA)

## Contents of `/`

|            |                                                         |
| :---       | :---                                                    |
| bin/       | Executables for running SCA analysis functions          |
| data/      | Input data (including those needed for the tutorials)   |
| notebooks/ | Example SCA notebooks                                   |
| Outputs/   | Output files from SCA analysis                          |
| pysca/     | Python code for SCA                                     |
| scripts/   | Utility scripts (e.g., getPfamDB.sh for Pfam database)  |

## Contents of `bin/`

|               |                                                                  |
| :---          | :---                                                             |
| annotateMSA   | Annotates alignments with phylogenetic/taxonomic information     |
| scaProcessMSA | Conducts some initial processing of the sequence alignment       |
| scaCore       | Runs the core SCA calculations                                   |
| scaSectorID   | Defines sectors given the results of the calculations in scaCore |

## Contents of `pysca/`

|             |                                                      |
| :---        | :---                                                 |
| scaTools.py | The SCA toolbox - functions for the SCA calculations |
| settings.py | Global configuration settings for the analysis       |

## Contents of `notebooks/`

Example notebooks demonstrating pySCA workflows will be available here. See [notebooks/README.md](notebooks/README.md) for details.
