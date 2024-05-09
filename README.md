# skin-digestion
Code necessary to reproduce analysis of human skin single cell RNA-seq datasets:
Clister, Terri, Rosalyn M. Fey, Zachary R. Garrison, Cristian D. Valenzuela, Anna Bar, Justin J. Leitenberger, and Rajan P. Kulkarni. 2024. "Optimization of Tissue Digestion Methods for Characterization of Photoaged Skin by Single Cell RNA Sequencing Reveals Preferential Enrichment of T Cell Subsets" Cells 13, no. 3: 266. https://doi.org/10.3390/cells13030266


Single cell RNA-seq and protein profiling datasets produced from non-affected human skin tissue collected from patients with skin cancer. Tissue was processed using one of two digestion methods.

File descriptions:
environment.yml
- yaml file used to create conda environment in which the analysis was run

environment_exported.yml
- yaml file containing all versions of all packages in conda environment used for these analyses
- exported from environment creating using the environment.yml file (above)

installRpackages.R
- R packages installed on the R terminal (not able to install using yml file)

integrateSCdata.R
- script used for analysis of fresh tissue samples (Figure 4)

fileList.RSEC_fresh.txt  
fileList.sampleTag_fresh.txt  
metadata_fresh_patientLevel.txt  
metadata_fresh.txt  
- input files for integrateSCdata.R

integrateSCdata_frozenVsFresh.R
- script used for analysis/comparison of fresh and frozen tissue samples

fileList.RSEC.txt  
fileList.SampleTag.txt  
metadata_patientLevel.txt  
metadata.txt  
- input files for integrateSCdata_frozenVsFresh.R
