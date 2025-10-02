This repository contains scripts to generate processed data and results in the following paper:

**Aberrant splicing in human cancer has a large-scale functional impact on transcription factors, Khalique Newaz, Olga Tsoy, and Jan Baumbach, under review.**

Two modes of usage are provided:

1. Reproducing the work in the paper
- Setup required: R version and packages need to be installed first. The code has been tested on R version 4.3.1.
- Install R version 4.3.1
- The following packages were used, which should be installed: "data.table version 1.17.8", "ggplot2 version 4.0.0", "ggrepel version 0.9.6", "dplyr version 1.1.4", "GenomicDataCommons version 1.26.0", "biomaRt  version 2.58.2", "seqinr version 4.2.36", "pheatmap version 1.0.13", "ggpp version 0.5.9", "TCGAbiolinks version 2.30.4"
- Run the script "runall_analysis.sh" as sh ./runall_analysis.sh on commandline terminal
- Two folders named "data" and "results_rep" with all processed data and results, respectively, will be created.
  
2. Redoing the analysis using newer data
- Here, some steps require manual 


