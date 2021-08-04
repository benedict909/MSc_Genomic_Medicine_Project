# PDAC Liquid Biopsy Bioinformatic Analysis Scripts

Scripts for the project. Run in numerical order to process .vcf files, generate inputs and then perform analysis. 

1. R script to load and combine .vcf files
2. a. shell script to install Annovar and download databases
2. b. shell script to annotate combined VCF from script 1 using Annovar
3. R script to filter variants in the annotated VCF
4. Shell scrip to convert .vcf files to .maf files for use with the R package `maftools`
5. R script to analyse combined VCF & MAF files for biomarkers
