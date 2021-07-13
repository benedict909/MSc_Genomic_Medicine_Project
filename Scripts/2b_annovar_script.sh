#!/bin/bash

# Liquid biopsies research project April 2021
# Run annovar

# Run ANNOVAR on VCF
# on "annovar_input" generated in the 1_load_concatanate_VCFs.R script
cp ~/../../mnt/c/Users/bened/Dropbox/Masters/9_Research_Project/Data/annovar_input.txt \
  ~/../../mnt/c/Users/bened/Dropbox/Masters/9_Research_Project/Data/annovar_input.avinput

todaysdate=$(date +%Y%m%d)

input_file=~/../../mnt/c/Users/bened/Dropbox/Masters/9_Research_Project/Data/annovar_input.avinput
output_file=~/../../mnt/c/Users/bened/Documents/Masters_research_project/annovar_workdir/annovar_output_all_variants_$todaysdate

# mkdir /mnt/c/Users/bened/Documents/Masters_research_project/annovar_workdir # make equivalent dir if you get "Error: cannot write LOG information to log file"

~/Annovar/annovar/table_annovar.pl $input_file ~/Annovar/annovar/humandb/ -buildver hg19 \
  -out $output_file -remove \
  -protocol refGene,cosmic94_coding,nci60,icgc28,clinvar_20210123,dbnsfp31a_interpro,avsnp150,gnomad211_genome,ALL.sites.2015_08 \
  -operation g,f,f,f,f,f,f,f,f -otherinfo -nastring . 
