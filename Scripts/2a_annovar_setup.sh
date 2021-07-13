#!/bin/bash

# Liquid biopsies research project April 2021
# Setup annovar
# N.B. both annovar and cosmic databases must be downloaded from the links below, as they are not opensource :(

mkdir ~/Annovar
cd ~/Annovar

# Then download annovar: https://www.openbioinformatics.org/annovar/annovar_download_form.php
tar -zxvf annovar.latest.tar.gz # N.B. this file must be downloaded from the Annovar website as it is not opensource

cd ~/Annovar/annovar 

# setup databases 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210123 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar nci60 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar icgc28 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
# other options (for full databases available see the annovar weebsite)
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/


# download COSMIC 94 & prepare database 
# get DL links from here: https://cancer.sanger.ac.uk/cosmic/download 
# Also check you are using the most recent version of COSMIC
myauth=$(echo "k20128353@kcl.ac.uk:BLUEeagle88!" | base64)

# download CosmicMutantExport
exportlink=$(curl -H "Authorization: Basic $myauth" \
  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v94/CosmicMutantExport.tsv.gz)

exportlink=$(echo $exportlink | sed 's/{"url"://g')
exportlink=$(echo $exportlink | sed 's/}//g')
exportlink=$(echo $exportlink | sed 's/"//g')

curl $exportlink --output CosmicMutantExport.tsv.gz

# download CosmicCodingMuts.vcf.gz
codinglink=$(curl -H "Authorization: Basic $myauth" \
  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v94/VCF/CosmicCodingMuts.vcf.gz)

codinglink=$(echo $codinglink | sed 's/{"url"://g')
codinglink=$(echo $codinglink | sed 's/}//g')
codinglink=$(echo $codinglink | sed 's/"//g')

curl $codinglink --output CosmicCodingMuts.vcf.gz
gunzip *.gz 

# create the database from the files downloaded from COSMIC
wget http://www.openbioinformatics.org/annovar/download/prepare_annovar_user.pl
chmod 775 prepare_annovar_user.pl
perl prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg19_cosmic94_coding.txt

mv hg19_cosmic94* ./humandb

rm CosmicMutantExport.tsv # remove to save disk space
rm CosmicCodingMuts.vcf