#! /bin/bash
set -eu

# Globals

pyscadir=pysca
scriptdir=scripts
datadir=data
outputdir=output

datarepo="https://github.com/ranganathanlab/pySCA-data"
version=6.0

# Download the data

cd ../

# In the event git is not installed, just directly download the data from
# GitHub using wget or curl (in order of preference). Also, check to see if tar
# is installed. If not, download the zipped archive.
if [ -x "$(command -v git)" ] && [ -d ".git/" ]; then
  git submodule init
  git submodule update --force
elif [ -x "$(command -v wget)" ]; then
  echo "git not installed --- trying wget"
  if [ -x "$(command -v tar)" ]; then
    wget -nc ${datarepo}/archive/v${version}.tar.gz
    tar xf v${version}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    wget -nc ${datarepo}/archive/v${version}.zip
    unzip v${version}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${datadir}
  mv -v pySCA-data-${version}/* ${datadir}/
  rmdir pySCA-data-${version}
elif [ -x "$(command -v curl)" ]; then
  echo "git not installed --- trying curl"
  if [ -x "$(command -v tar)" ]; then
    curl -L -O -C ${datarepo}/archive/v${version}.tar.gz
    tar xf v${version}.tar.gz
  elif [ -x "$(command -v unzip)" ]; then
    curl -L -O -C ${datarepo}/archive/v${version}.zip
    unzip v${version}.zip
  else
    echo "'unzip' or 'tar' (with gzip) is required for decompressing data."
    exit 3
  fi
  mkdir -p ${datadir}
  mv -v pySCA-data-${version}/* ${datadir}/
  rmdir pySCA-data-${version}
fi

# Generate the output files

mkdir -vp ${outputdir}

# The S1A serine proteases
echo "S1A serine protease Calculations:" | tee ${outputdir}/s1A_halabi.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/s1Ahalabi_1470_nosnakes.an -s 3TGI -c E \
  -t -n 2>&1 | tee -a ${outputdir}/s1A_halabi.log
python3 ${pyscadir}/scaCore.py ${outputdir}/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a ${outputdir}/s1A_halabi.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a ${outputdir}/s1A_halabi.log

# Beta-lactamase
echo "Beta-lactamase Calculations:" | tee ${outputdir}/PF13354.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/PF13354_full.an -s 1FQG -c A \
  -f 'Escherichia coli' -t -n 2>&1 | tee -a ${outputdir}/PF13354.log
python3 ${pyscadir}/scaCore.py ${outputdir}/PF13354_full.db 2>&1 | \
  tee -a ${outputdir}/PF13354.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/PF13354_full.db 2>&1 | \
  tee -a ${outputdir}/PF13354.log

# G-protein - this analysis is run with two alignments - the full PFAM
# alignment (PF00071_full) and the PFAM alignment filtered to remove several
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" | tee ${outputdir}/PF00071.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/PF00071_full.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n 2>&1 | tee -a ${outputdir}/PF00071.log
python3 ${pyscadir}/scaCore.py ${outputdir}/PF00071_full.db 2>&1 | \
  tee -a ${outputdir}/PF00071.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/PF00071_full.db 2>&1 | \
  tee -a ${outputdir}/PF00071.log

echo "G-protein calculations:" | tee ${outputdir}/PF00071_rd2.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/PF00071_rd2.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n 2>&1 | tee -a ${outputdir}/PF00071_rd2.log
python3 ${pyscadir}/scaCore.py ${outputdir}/PF00071_rd2.db 2>&1 | \
  tee -a ${outputdir}/PF00071_rd2.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/PF00071_rd2.db 2>&1 | \
  tee -a ${outputdir}/PF00071_rd2.log

# DHFR - this analysis is also run with two alignments for comparison -
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment
# (DHFR_PEPM3.an)
echo "DHFR Calculations:" | tee ${outputdir}/PF00186.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/PF00186_full.an -s 1RX2 -c A \
  -f 'Escherichia coli' -t -n 2>&1 | tee -a ${outputdir}/PF00186.log
python3 ${pyscadir}/scaCore.py ${outputdir}/PF00186_full.db 2>&1 | \
  tee -a ${outputdir}/PF00186.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/PF00186_full.db 2>&1 | \
  tee -a ${outputdir}/PF00186.log

echo "DHFR Calculations:" | tee ${outputdir}/DHFR_PEPM3.log
python3 ${pyscadir}/scaProcessMSA.py ${datadir}/DHFR_PEPM3.an -s 1RX2 -c A -t -n 2>&1 | \
  tee -a ${outputdir}/DHFR_PEPM3.log
python3 ${pyscadir}/scaCore.py ${outputdir}/DHFR_PEPM3.db 2>&1 | \
  tee -a ${outputdir}/DHFR_PEPM3.log
python3 ${pyscadir}/scaSectorID.py ${outputdir}/DHFR_PEPM3.db 2>&1 | \
  tee -a ${outputdir}/DHFR_PEPM3.log
