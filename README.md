# TCR-Cube
Bioinformatics investigation of dual-antigen recognition for the combination of TCR-T and Car-T cell therapies

## 🔬 Overview

TCR-T cells can recognise antigens expressed both on the cell surface and in intracellular compartments and TCR-T cell therapy has demonstrated encouraging potential for the treatment of solid tumours. However, cancer cells can resist the treatment by hiding the single target to prevent the recognition and survive, causing recurrence of cancers in patients. These challenges could be addressed by dual targeting strategies, the combination of pHLA and cell surface targeting T-cell therapies. Here we present a bioinformatics approach to identify dual target antigen pairs that address the tumour heterogeneity and antigen escape.

## 📁 Project structure
dual_target_pipeline/ \
├── README.md \
├── LICENSE \
├── .gitignore \
├── environment.yml               # (for conda env) or renv.lock if using R renv \
├── data/ \
│   ├── input/ \
│   │   ├── ct_antigens.csv \
│   │   ├── membrane_antigens.csv \
│   │   └── rna_seq_data.tsv      # (downloaded UCSC Toil TPM) \
│   └── processed/                # intermediate data \
├── scripts/ \
│   ├── 01_define_antigen_space.R \
│   ├── 02_extract_expression_data.R \
│   ├── 03_create_sketches.py \
│   ├── 04_train_antigen_pairs/ \
│   │   ├── identify_pairs.R \
│   │   └── run_all_indications.sh \
│   ├── 05_test_antigen_pairs/ \
│   │   ├── evaluate_pairs.R \
│   │   └── run_all_indications.sh \
│   └── 06_visualize_antigen_pair.R \
├── results/ \
│   ├── sketches/ \
│   ├── pair_identification/ \
│   ├── evaluation/ \
│   └── figures/ \
├── workflow/ \
│   └── run_pipeline.sh          # (master shell script) \
└── docs/ \
    └── method_overview.md 
