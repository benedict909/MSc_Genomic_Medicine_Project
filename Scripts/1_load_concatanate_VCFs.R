# Liquid Biopsy MSc project April 2021
# Load and combine VCF Files

library(vcfR)
library(tidyverse)

load2object <- function (filename) 
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error! the file ", filename, 
            " was not found :("))
  NULL
}

# Define directories 
raw_data_dir = "~/Masters_research_project/Raw_data" # directory containing all input .vcf files
per_sample_dir = "~/Masters_research_project/per_sample_vcf/"; if(!dir.exists(per_sample_dir)) dir.create(per_sample_dir)
vcf_dir = "~/Masters_research_project/VCF_files/"; if(!dir.exists(vcf_dir)) dir.create(vcf_dir)


# 1. Load each individual VCF and convert to dataframe ----

vcf_files = list.files(raw_data_dir)
# vcf_files = vcf_files[grepl("14",vcf_files)]
# filename = vcf_files[1]


for(filename in vcf_files){
  print(paste("loading and processing", filename))
  
  current_vcf = read.vcfR(file.path(raw_data_dir, filename))
  current_vcf =  vcfR2tidy(current_vcf, single_frame = T)
  current_vcf_df = current_vcf$dat
  rm(current_vcf) # to save memory
  
  if(grepl("CIRCB14",filename)){ # alternative method for this sample to include amount (e.g. "40ng")
    sample_notissue = "CIRCB14"
    if(grepl("PLASMA",filename)) tissue_name = "PLASMA"
    if(grepl("Blood",filename)) tissue_name = "Blood"
    amount =  gsub("_SmallVariants.genome.vcf","",unlist(strsplit(filename, tissue_name))[2])
    sample_name = paste(sample_notissue, tissue_name, amount, sep = "_")
  }else{
    filename_split = unlist(strsplit(filename, "_"))
    if(length(filename_split) == 2)  filename_split = unlist(strsplit(filename_split, split = "(?<=[0-9])", perl = T)) # fix for when there is no underscore between sample and tissue
    
    if(!is.na(as.numeric(filename_split[2]))){ # fix for when there are 2 numbers in the sample name 
      filename_split[1] = paste0(filename_split[1],filename_split[2])
      filename_split[2] = filename_split[3]
    }
    
    sample_notissue = filename_split[1]
    tissue_name = filename_split[2]
    if(grepl("CIRCB35_PLASMA",filename)) sample_notissue = "CIRCB35a"
    if(grepl("CIRCB35PLASMA",filename)) sample_notissue = "CIRCB35b"
    sample_name = paste(sample_notissue, tissue_name, sep = "_")
  }
  
  
  current_vcf_df = current_vcf_df %>% 
    mutate(Sample_notissue = sample_notissue,
           Tissue = tissue_name,
           Sample = sample_name) %>% 
    select(Sample, Sample_notissue, Tissue, everything()) %>% 
    as.data.frame()
  
  # save each individual samples VCF 
  save(current_vcf_df, file = paste0(per_sample_dir,sample_name,"_vcf.RData"))
}


# 2. load and combine each indvidual VCF RData ----

per_sample_files = list.files(per_sample_dir)
# filename = per_sample_files[1]

filter_PASS = c(TRUE, FALSE)[2] # select whether to filter VCF by FILTER == "PASS"

for(filename in per_sample_files){
  print(paste("Loading and adding", filename, "to the combined VCF"))
  
  current_vcf = load2object(file.path(per_sample_dir,filename))
  
  current_vcf_filt = filter(current_vcf,!is.na(ALT)) # keep variants only
  if(filter_PASS) current_vcf_filt = filter(current_vcf_filt, FILTER == "PASS") 
  
  rm(current_vcf) # to save memory
  
  if(filename == per_sample_files[1]){
    vcf_combined = current_vcf_filt
  }else{
    vcf_combined = rbind(vcf_combined, current_vcf_filt)
  }
}

# add extra info to VCF (N.B this must be changed manually) 
patient_1_samples = c("CIRCB1","CIRCB6","CIRCB19","CIRCB35","CIRCB35a", "CIRCB35b")
patient_2_samples = c("CIRCB4","CIRCB15","CIRCB40")
patient_3_samples = c("CIRCB8","CIRCB29")
patient_4_samples = c("CIRCB12","CIRCB23","CIRCB33")
patient_lung_samples = "CIRCB14"

timepoint_0_samples = c("CIRCB1","CIRCB4","CIRCB8","CIRCB12")
timepoint_1_samples = c("CIRCB6","CIRCB15","CIRCB23")
timepoint_2_samples = c("CIRCB19","CIRCB29","CIRCB33")
timepoint_3_samples = c("CIRCB35","CIRCB35a","CIRCB35b","CIRCB40")

# Add Patient/sample columns 
vcf_combined = vcf_combined %>% 
  mutate(Type = ifelse(nchar(ALT)==1 & nchar(REF)==1,"SNV",NA),
         Type = ifelse(nchar(ALT)==1 & nchar(REF)>1,"DEL",Type),
         Type = ifelse(nchar(ALT)>1 & nchar(REF)==1,"INS",Type),
         Patient = ifelse(  Sample_notissue %in% patient_1_samples, "Patient_1", NA),
         Patient = ifelse(  Sample_notissue %in% patient_2_samples, "Patient_2", Patient),
         Patient = ifelse(  Sample_notissue %in% patient_3_samples, "Patient_3", Patient),
         Patient = ifelse(  Sample_notissue %in% patient_4_samples, "Patient_4", Patient),
         Patient = ifelse(  Sample_notissue %in% patient_lung_samples, "Lung", Patient),
         Timepoint = ifelse(Sample_notissue %in% timepoint_0_samples, "T0", NA),
         Timepoint = ifelse(Sample_notissue %in% timepoint_1_samples, "T1", Timepoint),
         Timepoint = ifelse(Sample_notissue %in% timepoint_2_samples, "T2", Timepoint),
         Timepoint = ifelse(Sample_notissue %in% timepoint_3_samples, "T3", Timepoint),
         ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>% 
  select(Sample, Patient, Tissue, Timepoint, CHROM, POS, REF, ALT, Type, ID, everything()) %>% 
  select(-c(Indiv, Sample_notissue))
head(vcf_combined)

vcf_combined = vcf_combined %>% 
  mutate( ID = paste(CHROM, POS, REF, ALT, sep = "_"))

# fix for double base substitutions
vcf_dbs = vcf_combined %>% 
  filter(nchar(REF) == 2 & nchar(ALT) == 2) %>% 
  dplyr::slice(rep(row_number(), each = 2)) %>% 
  mutate(IDnew = paste0(Sample,ID)) %>% 
  mutate(duped = ifelse(duplicated(IDnew), T, F),
         REF = ifelse(!duped, substr(REF,1,1), REF),
         ALT = ifelse(!duped, substr(ALT,1,1), ALT),
         REF = ifelse(duped, substr(REF,2,2), REF),
         ALT = ifelse(duped, substr(ALT,2,2), ALT),
         POS = ifelse(duped, POS+1, POS),
         Type = "SNV") %>% 
  select(-c(duped,IDnew))

vcf_combined = vcf_combined %>% 
  filter(!(nchar(REF) == 2 & nchar(ALT) == 2)) %>% 
  rbind(vcf_dbs)

# fix for triple base subs 
vcf_trips = vcf_combined %>% 
  filter(is.na(Type)) %>% 
  dplyr::slice(rep(row_number(), each = 3)) %>% 
  mutate(trip_pos = rep(c(1:3),times = nrow(.)/3) ) %>% 
  mutate(REF = ifelse(trip_pos == 1, substr(REF,1,1), REF),
         ALT = ifelse(trip_pos == 1, substr(ALT,1,1), ALT),
         REF = ifelse(trip_pos == 2, substr(REF,2,2), REF),
         ALT = ifelse(trip_pos == 2, substr(ALT,2,2), ALT),
         POS = ifelse(trip_pos == 2, POS+1, POS),
         REF = ifelse(trip_pos == 3, substr(REF,3,3), REF),
         ALT = ifelse(trip_pos == 3, substr(ALT,3,3), ALT),
         POS = ifelse(trip_pos == 3, POS+2, POS),
         Type = "SNV") %>% 
  select(-trip_pos)

vcf_combined = vcf_combined %>% 
  filter(!is.na(Type)) %>% 
  rbind(vcf_trips) %>% 
  mutate(ID = paste(CHROM, POS, REF, ALT, sep = "_"))


# save combined VCFs 
if(filter_PASS)  combined_filename = "VCF_allsamples_filtered_allcolumns.RData"
if(!filter_PASS) combined_filename = "VCF_allsamples_unfiltered_allcolumns.RData"

save(vcf_combined, file = file.path(vcf_dir,combined_filename))
save(vcf_combined, file = file.path("~/../Dropbox/Masters/9_Research_Project/Data",combined_filename))

# load(file = file.path(vcf_dir,combined_filename))

# create annovar input 
annovar_input = vcf_combined %>% 
  filter(Type =="SNV") %>% 
  select(CHROM, Start= POS, End = POS, REF, ALT) %>% 
  distinct()

write_delim(x = annovar_input, file = "~/../Dropbox/Masters/9_Research_Project/Data/annovar_input.txt",
            delim = "\t", col_names = F)
