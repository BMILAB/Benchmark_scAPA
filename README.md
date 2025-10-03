# Benchmarking computational methods for identifying and quantifying polyadenylation sites from single-cell RNA-seq data
This repository contains R scripts that were used to generate results and plots for the paper titled "**Benchmarking computational methods for identifying and quantifying polyadenylation sites from single-cell RNA-seq data**".  

**1. DataProcessing/**: This directory contains scripts for data pre-processing of the benchmark study.
  
  **_functions.R_**: Common functions called by other scripts.
  
  **_readPACs.R_**: This script is used to  reads the outputs from APA analysis methods and stores them as a *PACdataset* object from [movAPA](https://github.com/BMILAB/movAPA/) for filtering, statistics, and other operations.
  
  **_sata_pa.R_**: This script is used to integrate *PACdataset* objects from different methods for the same sample and store them as a list for subsequent unified data loading and calculation. Meanwhile, it preprocesses the reference pA dataset.
  
  
**2. Benchmark/**: This directory contains scripts to evaluate the performance of the ten APA analysis methods under various aspects.

  **_identification_benchmark.R_**: This script is used to compare the pAs identified by different methods and generate the corresponding figures for the manuscript.
  
  **_uniquePA_benchmark.R_**: This script is used to compare the consistency of pA identification results across different methods and generate the corresponding figures for the manuscript.
  
  **_quantification_benchmark.R_**: This script is used to evaluate the performance of different methods for pA quantification and generate the corresponding figures for the manuscript.
  
  **_DEAPA_benchmark.R_**: This script is used to evaluate the performance of different methods for detecting DEAPA genes and generate the corresponding figures for the manuscript.
  
  **_summary_benchmark.R_**: This script is used to compare the time and memory consumption of different methods, summarize evaluation results, and generate corresponding figures for the manuscript.
  
