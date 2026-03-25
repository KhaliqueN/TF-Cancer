#!/usr/bin/sh

wget https://zenodo.org/records/17250324/files/data_for_study_reproduction.zip
unzip data_for_study_reproduction.zip

Rscript TF_DBD_annotations.r
Rscript TF_DBD_splicing_patterns_annotations.r
Rscript TF_ED_splicing_patterns_annotations.r
sh ./run_reproduce.sh
