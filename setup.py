#! /usr/bin/env python3
""" Installation and setup for pySCA"""

from setuptools import setup

setup(
    name="pySCA",
    version="6.1",
    author="Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds",
    maintainer="Ansel George",
    packages=["pysca"],
    package_data={"pysca": ["settings.py"]},
    description="Python 3 implementation of Statistical Coupling Analysis (SCA)",
    url="https://ranganathanlab.gitlab.io/pySCA",
    download_url="https://github.com/ranganathanlab/pySCA",
    long_description=open("README.md", "r").read(),
    install_requires=[
        "biopython",
        "numpy",
        "scipy",
        "argparse",
        "wheel",
        "matplotlib",
    ],
    extras_require={
        "notebooks": [
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",  # Required for Jupyter kernel registration
            "ipywidgets>=7.0.0",
            "plotly>=5.0.0",  # For interactive visualizations (primary)
            "bokeh>=2.0.0",  # Alternative visualization library
        ],
    },
    scripts=[
        "bin/alnChangeDelim",
        "bin/alnFilterSeqSize",
        "bin/alnParseID",
        "bin/annotateMSA",
        "bin/scaProcessMSA",
        "bin/sca-process-msa",  # Python 3 wrapper
        "bin/alnConvertGI",
        "bin/alnReplaceHeaders",
        "bin/scaCore",
        "bin/sca-core",  # Python 3 wrapper
        "bin/scaSectorID",
        "bin/sca-sectorid",  # Python 3 wrapper
    ],
)
