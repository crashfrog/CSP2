#!/usr/bin/env bash
set -euo pipefail

BASE=/home/Justin.Payne

# ---- 1. MUMmer baseline ----
mkdir -p $BASE/snpdiffs_cmp/mummer
cd $BASE/snpdiffs_cmp/mummer
dnadiff -p Q1_vs_Q2 $Q1 $Q2 2>&1 | grep -E "^[0-9]"
echo "===MUMMER_SNP==="
show-snps -T Q1_vs_Q2.delta 2>&1
echo "===MUMMER_COORDS==="
show-coords -T -c Q1_vs_Q2.delta 2>&1

# ---- 2. Phraya ----
mkdir -p $BASE/snpdiffs_cmp/phraya
cd $BASE/snpdiffs_cmp/phraya

# phraya plan: query1 is the reference, query2 is the input
$($PHRAYA plan --reference $Q1 --inputs $Q2 --output plan.phrayaplan 2>&1)
echo "===PLAN_TASKS==="
$PHRAYA plan-tasks plan.phrayaplan 2>&1
echo "===ALIGN==="
$PHRAYA align --reference $Q1 plan.phrayaplan --output aln_out --strategy balanced 2>&1
echo "===PHRAYA_FILE==="
ls aln_out/ 2>&1
echo "===FILTER_SNPDIFFS==="
$PHRAYA filter aln_out/ref.phraya --format snpdiffs --output pair.snpdiffs --reference-fasta $Q1 --query-fasta $Q2 --reference-id "query1" 2>&1
echo "===SNPDIFFS_OUTPUT==="
cat pair.snpdiffs
