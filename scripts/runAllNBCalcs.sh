#! /bin/bash
set -eu

cd ../

git submodule init
git submodule update --force

pyscadir=pysca
scriptdir=scripts
datadir=data
outputdir=output

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
