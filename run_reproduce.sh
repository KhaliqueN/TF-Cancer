#!/usr/bin/sh

wget https://zenodo.org/records/17250324/files/data_for_study_reproduction.zip
unzip data_for_study_reproduction.zip

Rscript TCGASpliceSeq_data.r
Rscript TF_splicing_patterns.r
Rscript TF_DBD_splicing_patterns_analysis.r
Rscript TF_ED_splicing_patterns_analysis.r
Rscript TF_dependency.r
Rscript TF_ATACseq.r
Rscript TF_Lambourne_annotations.r
Rscript map_hg38.r
