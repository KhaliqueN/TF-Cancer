This repository contains scripts to generate processed data and results in the following paper:

**Aberrant splicing in human cancer has a large-scale functional impact on transcription factors, Khalique Newaz, Olga Tsoy, and Jan Baumbach, under review.**

The code has been tested wtih R version 4.3.1 and the following package versions, which should be installed before runnign the scirpts: "data.table version 1.17.8", "ggplot2 version 4.0.0", "ggrepel version 0.9.6", "dplyr version 1.1.4", "GenomicDataCommons version 1.26.0", "biomaRt  version 2.58.2", "seqinr version 4.2.36", "pheatmap version 1.0.13", "ggpp version 0.5.9", "TCGAbiolinks version 2.30.4", "Biostrings version 2.70.3"

The repository can be used to two ways:

1. Reproducing the work in the paper
- Run the script "run_reproduce.sh" as sh ./run_reproduce.sh on commandline terminal
- Two folders named "data" and "results_rep" with all processed data and results, respectively, will be created.
  
2. Producing the results of the work using newer UniProt versions
- Run the script "run_with_current_UniProt.sh" as sh ./run_with_current_UniProt.sh on commandline terminal
- Two folders named "data" and "results_rep" with all processed data and results, respectively, will be created.


