#!/bin/bash

# install vcf2maf tool from https://github.com/mskcc/vcf2maf

# # N.B you will need samtools installed
# # check https://repo.anaconda.com to find the latest version 
# # cd ~/
# wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh # Install Anaconda 
# chmod +x ./Anaconda3-2021.05-Linux-x86_64.sh
# bash ./Anaconda3-2021.05-Linux-x86_64.sh
# source ~/.bashrc
# conda install -c samtools # Install required packages with Anaconda 
# # conda install --override-channels -c main -c bioconda samtools # alternative if conda install not working

# export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
# curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*

# # install VEP
# # conda create --name myenv2 # if installation doesn't work, try this
# # conda activate myenv2
# conda install -qy -c conda-forge -c bioconda -c defaults ensembl-vep==102.0 htslib==1.10.2 bcftools==1.10.2 samtools==1.10 ucsc-liftover==377

# # install ref genome for VEP
# mkdir -p $HOME/.vep/homo_sapiens/102_GRCh37/
# rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh37.tar.gz $HOME/.vep/
# tar -zxf $HOME/.vep/homo_sapiens_vep_102_GRCh37.tar.gz -C $HOME/.vep/
# rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/grch37/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz $HOME/.vep/homo_sapiens/102_GRCh37/
# gzip -d $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
# bgzip -i $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa
# samtools faidx $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz

#conda init bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate myenv2

# 1. run conversion on somatic VCFs
maf_dir=~/../../mnt/c/Users/bened/Documents/Masters_research_project/MAF_input
maf_resdir=~/../../mnt/c/Users/bened/Documents/Masters_research_project/MAF_files

mkdir $maf_dir/MAF_input_vcf
mkdir $maf_resdir

# vcf colnames string
colnames_str="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL"
vcf_suffix=_vcf_filtered.txt


input_files=$(ls $maf_dir/MAF_input_txt)

for file in $input_files; do
	echo $file 
	samp_name=$(echo $file | sed 's/_vcf_filtered.txt//g')
	
	# create input file 
	maf_input_file=$maf_dir/MAF_input_vcf/$samp_name.maf_input.vcf
	touch $maf_input_file
	cat ~/../../mnt/c/Users/bened/Documents/Masters_research_project/Raw_data/CIRCB1_PLASMA_SmallVariants.genome.vcf | head -n 50 > $maf_input_file
	echo -e $colnames_str >> $maf_input_file
	cat $maf_dir/MAF_input_txt/$samp_name$vcf_suffix >> $maf_input_file

	# name of output filename
	output_filepath=$maf_resdir/$samp_name.maf

	# run VCF to MAF conversion
	perl ~/mskcc-vcf2maf-754d68a/vcf2maf.pl \
	  --input-vcf $maf_input_file \
	  --output-maf $output_filepath \
	  --ref-fasta ~/../../mnt/c/Users/bened/Documents/Genomes/hg19.fa \
	  --vep-path ~/anaconda3/envs/myenv2/bin \
	  --tumor-id $samp_name --normal-id CIRCB1_Blood \
	  --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --vep-overwrite
done

# 2. run conversion on germline VCFs

mkdir $maf_dir/MAF_input_blood_vcf

input_files=$(ls $maf_dir/MAF_input_blood_txt)

for file in $input_files; do
	echo $file 
	samp_name=$(echo $file | sed 's/_vcf_filtered.txt//g')
	
	# create input file 
	maf_input_file=$maf_dir/MAF_input_blood_vcf/$samp_name.maf_input.vcf
	touch $maf_input_file
	cat ~/../../mnt/c/Users/bened/Documents/Masters_research_project/Raw_data/CIRCB1_PLASMA_SmallVariants.genome.vcf | \
		head -n 50 \
		> $maf_input_file
	echo $colnames_str >> $maf_input_file
	vcf_suffix=_vcf_filtered.txt
	cat $maf_dir/MAF_input_blood_txt/$samp_name$vcf_suffix >> $maf_input_file

	# name of output filename
	output_filename=$samp_name.maf
	output_filepath=~/../../mnt/c/Users/bened/Documents/Masters_research_project/MAF_files/$samp_name.maf

	# run VCF to MAF conversion
	perl ~/mskcc-vcf2maf-754d68a/vcf2maf.pl \
	  --input-vcf $maf_input_file \
	  --output-maf $output_filepath \
	  --ref-fasta ~/../../mnt/c/Users/bened/Documents/Genomes/hg19.fa \
	  --vep-path ~/anaconda3/envs/myenv2/bin \
	  --tumor-id $samp_name --normal-id CIRCB1_Blood \
	  --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --vep-overwrite
done
