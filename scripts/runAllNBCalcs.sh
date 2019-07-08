#! /bin/bash
set -eu

cd ../

git submodule init
git submodule update --force

mkdir -vp output

# The S1A serine proteases
echo "S1A serine protease Calculations:" | tee output/s1A_halabi.log
python3 pysca/scaProcessMSA.py data/s1Ahalabi_1470_nosnakes.an -s 3TGI -c E \
  -t -n 2>&1 | tee -a output/s1A_halabi.log
python3 pysca/scaCore.py output/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a output/s1A_halabi.log
python3 pysca/scaSectorID.py output/s1Ahalabi_1470_nosnakes.db 2>&1 | \
  tee -a output/s1A_halabi.log

# Beta-lactamase
echo "Beta-lactamase Calculations:" | tee output/PF13354.log
python3 pysca/scaProcessMSA.py data/PF13354_full.an -s 1FQG -c A \
  -f 'Escherichia coli' -t -n 2>&1 | tee -a output/PF13354.log
python3 pysca/scaCore.py output/PF13354_full.db 2>&1 | \
  tee -a output/PF13354.log
python3 pysca/scaSectorID.py output/PF13354_full.db 2>&1 | \
  tee -a output/PF13354.log

# G-protein - this analysis is run with two alignments - the full PFAM
# alignment (PF00071_full) and the PFAM alignment filtered to remove several
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" | tee output/PF00071.log
python3 pysca/scaProcessMSA.py data/PF00071_full.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n 2>&1 | tee -a output/PF00071.log
python3 pysca/scaCore.py output/PF00071_full.db 2>&1 | \
  tee -a output/PF00071.log
python3 pysca/scaSectorID.py output/PF00071_full.db 2>&1 | \
  tee -a output/PF00071.log

echo "G-protein calculations:" | tee output/PF00071_rd2.log
python3 pysca/scaProcessMSA.py data/PF00071_rd2.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n 2>&1 | tee -a output/PF00071_rd2.log
python3 pysca/scaCore.py output/PF00071_rd2.db 2>&1 | \
  tee -a output/PF00071_rd2.log
python3 pysca/scaSectorID.py output/PF00071_rd2.db 2>&1 | \
  tee -a output/PF00071_rd2.log

# DHFR - this analysis is also run with two alignments for comparison -
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment
# (DHFR_PEPM3.an)
echo "DHFR Calculations:" | tee output/PF00186.log
python3 pysca/scaProcessMSA.py data/PF00186_full.an -s 1RX2 -c A \
  -f 'Escherichia coli' -t -n 2>&1 | tee -a output/PF00186.log
python3 pysca/scaCore.py output/PF00186_full.db 2>&1 | \
  tee -a output/PF00186.log
python3 pysca/scaSectorID.py output/PF00186_full.db 2>&1 | \
  tee -a output/PF00186.log

echo "DHFR Calculations:" | tee output/DHFR_PEPM3.log
python3 pysca/scaProcessMSA.py data/DHFR_PEPM3.an -s 1RX2 -c A -t -n 2>&1 | \
  tee -a output/DHFR_PEPM3.log
python3 pysca/scaCore.py output/DHFR_PEPM3.db 2>&1 | \
  tee -a output/DHFR_PEPM3.log
python3 pysca/scaSectorID.py output/DHFR_PEPM3.db 2>&1 | \
  tee -a output/DHFR_PEPM3.log
