#!/bin/bash
# Master script to run the full dual-target pipeline

set -e

# Step 1: Extract RNA-seq data
Rscript scripts/01_extract_expression_data.R

# Step 2: Create geometric sketch
python3 scripts/02_create_sketches.py

# Step 3: Train antigen pair classifier
seq 1 33 | xargs -P 8 -n 1 bash scripts/03_train_antigen_pairs/run_all_indications.sh

# Step 4: Evaluate on test set
seq 1 33 | xargs -P 8 -n 1 bash scripts/04_test_antigen_pairs/run_all_indications.sh
