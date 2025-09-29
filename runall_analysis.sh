#!/usr/bin/sh

Rscript preprocess_data.r
Rscript TF_DBD_annotations.r
Rscript TCGASpliceSeq_data.r
Rscript TF_splicing_patterns.r
Rscript TF_DBD_splicing_patterns_annotations.r
Rscript TF_DBD_splicing_patterns_analysis.r
Rscript TF_ED_splicing_patterns_annotations.r
Rscript TF_ED_splicing_patterns_analysis.r
Rscript TF_dependency.r
Rscript TF_ATACseq.r
