# TCR-Cube
Bioinformatics investigation of dual-antigen recognition for the combination of TCR-T and Car-T cell therapies

## 🔬 Overview

1. **Define antigen search space** from CT and membrane gene lists.
2. **Extract RNA-seq TPM values** from UCSC Toil-processed TCGA/GTEx data.
3. **Subsample with geometric sketching** to balance tissue distribution.
4. **Identify antigen pairs** with separation of tumor vs normal (DB index).
5. **Evaluate performance** on test dataset.
6. **Visualize top pairs** for presentation and validation.

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
