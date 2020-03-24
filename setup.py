#! /usr/bin/env python3
""" Installation and setup for pySCA"""

from setuptools import setup

setup(
    name="pySCA",
    version="6.1",
    author="Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds, and Ansel George",
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
        "jupyterlab",
    ],
    scripts=[
        "bin/alnChangeDelim",
        "bin/alnFilterSeqSize",
        "bin/alnParseID",
        "bin/annotateMSA",
        "bin/scaProcessMSA",
        "bin/alnConvertGI",
        "bin/alnReplaceHeaders",
        "bin/scaCore",
        "bin/scaSectorID",
    ],
)
