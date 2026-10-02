#!/usr/bin/env bash
set -euo pipefail

export PATH=/nfs/software/modules/mummer/4.0.0/bin:$PATH
export HOME=/home/Justin.Payne
export CARGO_HOME=$HOME/.cargo

PHRAYA=$HOME/phraya_src/target/release/phraya
Q1=$HOME/phraya_src/query1.fasta
Q2=$HOME/phraya_src/query2.fasta

echo "===Q1 header==="; head -1 $Q1
echo "===Q2 header==="; head -1 $Q2

mkdir -p $HOME/snph2
cd $HOME/snph2
rm -rf plan.phrayaplan aln_out pair.snpdiffs

echo "===PLAN==="
$PHRAYA plan --reference $Q1 --inputs $Q2 --output plan.phrayaplan 2>&1
echo "===ALIGN==="
$PHRAYA align --reference $Q1 plan.phrayaplan --output aln_out --strategy balanced 2>&1
echo "===FILES==="
ls aln_out/

# Find the .phraya file (named after reference contig)
REFPHRAYA=$(ls aln_out/*.phraya 2>/dev/null | head -1)
echo "===REF_PHARYA===$REFPHRAYA"

echo "===FILTER==="
$PHRAYA filter "$REFPHRAYA" --format snpdiffs --output pair.snpdiffs --reference-fasta $Q1 --query-fasta $Q2 --reference-id "query1" --query-id "query2" 2>&1

echo "===OUTPUT==="
cat pair.snpdiffs
echo "===END==="
