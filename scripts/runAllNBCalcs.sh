#! /bin/bash
set -eu

cd ../

git submodule init
git submodule update --force

mkdir -vp output

# The S1A serine proteases
echo "S1A serine protease Calculations:" > output/s1A_halabi.log 2>&1
python3 pysca/scaProcessMSA.py data/s1Ahalabi_1470_nosnakes.an -s 3TGI -c E -t -n \
  >> output/s1A_halabi.log 2>&1
python3 pysca/scaCore.py output/s1Ahalabi_1470_nosnakes.db \
  >> output/s1A_halabi.log 2>&1
python3 pysca/scaSectorID.py output/s1Ahalabi_1470_nosnakes.db \
  >> output/s1A_halabi.log 2>&1

# Beta-lactamase
echo "Beta-lactamase Calculations:" > output/PF13354.log
python3 pysca/scaProcessMSA.py data/PF13354_full.an -s 1FQG -c A \
  -f 'Escherichia coli' -t -n >> output/PF13354.log 2>&1
python3 pysca/scaCore.py output/PF13354_full.db >> output/PF13354.log 2>&1
python3 pysca/scaSectorID.py output/PF13354_full.db >> output/PF13354.log 2>&1

# G-protein - this analysis is run with two alignments - the full PFAM
# alignment (PF00071_full) and the PFAM alignment filtered to remove several
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" > output/PF00071.log
python3 pysca/scaProcessMSA.py data/PF00071_full.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n >> output/PF00071.log 2>&1
python3 pysca/scaCore.py output/PF00071_full.db >> output/PF00071.log 2>&1
python3 pysca/scaSectorID.py output/PF00071_full.db >> output/PF00071.log 2>&1

echo "G-protein calculations:" > output/PF00071_rd2.log
python3 pysca/scaProcessMSA.py data/PF00071_rd2.an -s 5P21 -c A \
  -f 'Homo sapiens' -t -n >> output/PF00071_rd2.log 2>&1
python3 pysca/scaCore.py output/PF00071_rd2.db >> output/PF00071_rd2.log 2>&1
python3 pysca/scaSectorID.py output/PF00071_rd2.db >> output/PF00071_rd2.log 2>&1

# DHFR - this analysis is also run with two alignments for comparison -
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment
# (DHFR_PEPM3.an)
echo "DHFR Calculations:" > output/PF00186.log
python3 pysca/scaProcessMSA.py data/PF00186_full.an -s 1RX2 -c A \
  -f 'Escherichia coli' -t -n >> output/PF00186.log 2>&1
python3 pysca/scaCore.py output/PF00186_full.db >> output/PF00186.log 2>&1
python3 pysca/scaSectorID.py output/PF00186_full.db >> output/PF00186.log 2>&1

echo "DHFR Calculations:" > output/DHFR_PEPM3.log
python3 pysca/scaProcessMSA.py data/DHFR_PEPM3.an -s 1RX2 -c A -t -n \
  >> output/DHFR_PEPM3.log 2>&1
python3 pysca/scaCore.py output/DHFR_PEPM3.db >> output/DHFR_PEPM3.log 2>&1
python3 pysca/scaSectorID.py output/DHFR_PEPM3.db >> output/DHFR_PEPM3.log 2>&1
