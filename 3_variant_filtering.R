# Script 3 - Variant filtering on annotated VCF & prepare MAF input

library(tidyverse)
library(vcfR)

load2object <- function (filename) {
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error! the file ", filename, 
            " was not found :("))
  NULL
}


# Load vcf_combined
vcf_dir = "~/Masters_research_project/VCF_files"
load(file = file.path(vcf_dir,"VCF_allsamples_unfiltered_allcolumns.RData"))



# 1. load annovar annots ----
annovar_load = read.delim2(paste0("~/Masters_research_project/annovar_workdir/annovar_output_all_variants_",
                                  #gsub("-","",Sys.Date()), ".hg19_multianno.txt"),na.strings = ".")
                                  "20210616.hg19_multianno.txt"),na.strings = ".")
save(annovar_load, file = file.path(vcf_dir,"annovar_VCF_allsamples_unfiltered.RData"))

keepcols =  colnames(annovar_load)[c(6:22,39)] # keep these annotation columns from the annovar output

annovar_annot = annovar_load %>% 
  distinct() %>% 
  mutate(ID = paste(Chr, Start, Ref, Alt, sep = "_")) %>% 
  select(ID, all_of(keepcols)) %>% 
  rename(gnomAD_AF = AF,  "1000G" = ALL.sites.2015_08)

# combine annovar output with vcf 
vcf_annotated = vcf_combined %>% 
  left_join(annovar_annot) %>% 
  separate(gt_AD, c("ref_DP","alt_DP")) %>%
  mutate(Sample_ID = paste(Patient, ID, sep = "_"), 
         ref_DP = as.numeric(ref_DP),
         alt_DP = as.numeric(alt_DP),
         gt_VF = as.numeric(gt_VF),
         Sample = ifelse(Sample == "CIRCB14_Blood_","CIRCB14_Blood",Sample))

# save vcf_annotated
save(vcf_annotated, file = file.path(vcf_dir,"VCF_allsamples_unfiltered_allcolumns_annovar.RData"))
save(vcf_annotated, file = file.path("~/../Dropbox/Masters/9_Research_Project/Data/VCF_allsamples_filtered_allcolumns_annovar.RData"))


# 2. variant filtering ----

# germline variants
vcf_blood = vcf_annotated %>% 
  filter(Tissue == "Blood" & FILTER == "PASS") 

# PDAC somatic variants = variants not found in blood
blood_IDS = vcf_blood$Sample_ID

vcf_somatic = vcf_annotated %>% 
  filter(FILTER == "PASS" & Tissue != "Blood" & !Sample_ID %in% blood_IDS)

vcf_somatic %>% # no. of filtered somatic variants per sample 
  count(Sample, Timepoint, Patient) %>% 
  arrange(Patient, Timepoint, n)

# add back in LowSupport variants seen in other samples of the same patient
somatic_IDs = vcf_somatic$Sample_ID 

vcf_annotated %>% 
  filter(FILTER == "LowSupport" & Tissue != "Blood" & 
           Sample_ID %in% somatic_IDs,
         alt_DP < 30) %>% 
  ggplot(aes(x = alt_DP)) +
    geom_bar(stat = "count") + 
    theme_bw() +  scale_x_continuous(minor_breaks = c(seq(0,30,1)), breaks = c(seq(0,30,5))) + 
    labs(y = "Count", x = "Alteration depth", title = "Alteration depth of `LowSupport` variants in plasma/urine")

vcf_annotated %>% 
  filter(FILTER == "PASS" & Tissue != "Blood" & 
           Sample_ID %in% somatic_IDs,
         alt_DP < 30) %>% 
  ggplot(aes(x = alt_DP)) +
    geom_bar(stat = "count") + 
    theme_bw() + scale_x_continuous(minor_breaks = c(seq(0,30,1)), breaks = c(seq(0,30,5))) + 
    labs(y = "Count", x = "Alteration depth", title = "Alteration depth of `PASS` variants in plasma/urine")

vcf_restore = vcf_annotated %>% 
  filter(FILTER == "LowSupport" & Tissue != "Blood" & 
           Sample_ID %in% somatic_IDs & alt_DP >= 3) # keep lowsupport variants with at least 3 ALT observations

vcf_restore %>% 
  ggplot(aes(x = QUAL)) +
    geom_bar() + # a Pred quality score of 20 = 99% base call accuracy 
    theme_bw() + 
    labs(title = "Sequencing quality of `LowSupport` variants with alt_DP >= 3")

# combine lowsupport variants with ALT_DP >= 3 with PASS variants 
vcf_somatic = rbind(vcf_somatic, vcf_restore)
count(vcf_somatic,Sample,Timepoint, Patient) %>% # new no. of somatic variants per sample 
  arrange(Patient, Timepoint, n)

# print some stats
table(vcf_somatic$Func.refGene) # no. of exonic variant in all samples
table(vcf_somatic$ExonicFunc.refGene) # no. of non-synonymous exonic SNVs in all samples 
table(vcf_somatic$Type) # variant type

# a greater proportion of exonic variants are lost to filtering than intronic 
sort(table(vcf_annotated$Func.refGene) / nrow(vcf_annotated))
sort(table(vcf_somatic$Func.refGene) / nrow(vcf_somatic))

# how many germline variants are seen in LBs 
vcf_lb = vcf_annotated %>% 
  filter(FILTER == "PASS" & Tissue != "Blood")
table(blood_IDS %in% vcf_lb$Sample_ID) # only 68 gl vars not seen in the liquid biopsies


save(vcf_somatic, file = file.path(vcf_dir,"vcf_LBsamples_filtered_allcolumns.RData"))
save(vcf_somatic, file = file.path("~/../Dropbox/Masters/9_Research_Project/Data/vcf_LBsamples_filtered_allcolumns_pluslung.RData"))


# 3. MAF input generate ----
# The rest of this script is only necessary if you want to perform analyses with maftools (e.g. make oncoplots)
# The original VCF files must be used as MAF input requires tumour/normal genotype info which is only in the original VCF, not vcf_annotated

# Define directories for MAF
raw_data_dir = "~/../Documents/Masters_research_project/Raw_data" # containing initial .vcf  files
per_sample_dir = "~/Masters_research_project/per_sample_vcf/" # containing _vcf.RData files
maf_input_dir = "~/Masters_research_project/MAF_input/MAF_input_txt"; if(!file.exists(maf_input_dir)) dir.create(maf_input_dir)
maf_input_dir_blood = "~/Masters_research_project/MAF_input/MAF_input_blood_txt"; if(!file.exists(maf_input_dir_blood)) dir.create(maf_input_dir_blood)

load(file.path(vcf_dir,"vcf_LBsamples_filtered_allcolumns.RData"))

# somatic MAF generate
somatic_IDs = vcf_somatic$Sample_ID
vcf_files = list.files(raw_data_dir)#; vcf_files = vcf_files[!grepl("vep",vcf_files)]

for(patient in sort(unique(vcf_annotated$Patient))){
  print(patient)
  current_samps = unique(vcf_annotated$Sample[vcf_annotated$Patient == patient])
  current_samps = c(current_samps[grepl("Blood", current_samps)], # reorder to put blood sample 1st
                    current_samps[!grepl("Blood", current_samps)])
  
  for(samp in current_samps){
    print(samp)
    
    samp_alt = gsub("_","",samp)
    mypattern = paste0(samp,"|",samp_alt)
    
    if(mypattern == "CIRCB35a_PLASMA|CIRCB35aPLASMA") mypattern = "CIRCB35_PLASMA";samp =  "CIRCB35a_PLASMA"
    if(mypattern == "CIRCB35b_PLASMA|CIRCB35bPLASMA") mypattern = "CIRCB35PLASMA"; samp = "CIRCB35b_PLASMA"
    
    current_load = NA
    current_load = vcf_files[grepl(mypattern, vcf_files)]
    
    for(i in 1:length(current_load)){
      if(is.na(current_load[i])) stop(paste(samp, "was not found in vcf_files :("))
      current_vcf = read.vcfR(file.path(raw_data_dir, current_load[i]))
      
      if(grepl("Blood", samp, ignore.case = T)){
        normal_gt_df = current_vcf@gt
        colnames(normal_gt_df)[2] = "NORMAL"
        
        normal_info = getINFO(current_vcf) %>%
          as.data.frame() %>% 
          rename(INFO = ".")
        
        normal_df = getFIX(current_vcf) %>%
          as.data.frame() %>% 
          mutate(ID = paste(patient, CHROM, POS, substr(REF,1,1), sep = "_")) %>% 
          cbind(normal_info, normal_gt_df) %>% 
          select(ID, NORMAL)
        
      }else{
        tumour_gt_df = current_vcf@gt
        colnames(tumour_gt_df)[2] = "TUMOR"
        
        tumour_info = getINFO(current_vcf) %>%
          as.data.frame() %>% 
          rename(INFO = ".")
        
        tumour_df = getFIX(current_vcf) %>%
          as.data.frame() %>% 
          mutate(ID = paste(patient, CHROM, POS, substr(REF,1,1), sep = "_")) %>% 
          cbind(tumour_info, tumour_gt_df) %>% 
          left_join(normal_df) %>% 
          filter(!is.na(ALT))
        
        tumour_df_filt = tumour_df %>% 
          mutate(ID = paste(patient, CHROM, POS, REF, ALT, sep = "_")) %>% 
          filter(ID %in% somatic_IDs)
        
        write_delim(tumour_df_filt, 
                    file = file.path(maf_input_dir,paste0(samp,"_vcf_filtered.txt")),
                    delim = "\t",col_names = F)
      }
    }
  }
}


# germline MAF generate
# load(file.path(vcf_dir,"VCF_allsamples_unfiltered_allcolumns_annovar.RData"))
gl_IDs = vcf_annotated %>% 
  filter(Tissue == "Blood", FILTER == "PASS", !is.na(ALT)) %>% 
  pull(Sample_ID)

vcf_files = list.files(raw_data_dir)#; vcf_files = vcf_files[!grepl("vep",vcf_files)]

for(patient in sort(unique(vcf_annotated$Patient))){
  print(patient)

  samp = unique(vcf_annotated$Sample[vcf_annotated$Patient == patient & vcf_annotated$Tissue == "Blood"]) 
  samp_alt = gsub("_","",samp)
  mypattern = paste0(samp,"|",samp_alt)
  
  print(samp)
  
  current_load = NA
  current_load = vcf_files[grepl(mypattern, vcf_files)]
  if(is.na(current_load)) stop(paste(samp, "was not found in vcf_files :("))
  
  current_vcf = read.vcfR(file.path(raw_data_dir, current_load))
  
  normal_gt_df = current_vcf@gt
  colnames(normal_gt_df)[2] = "NORMAL"
  
  normal_info = getINFO(current_vcf) %>%
    as.data.frame() %>% 
    rename(INFO = ".")
  
  normal_df = getFIX(current_vcf) %>%
    as.data.frame() %>% 
    mutate(ID = paste(samp, CHROM, POS, REF, sep = "_")) %>% 
    cbind(normal_info, normal_gt_df) %>% 
    select(ID, NORMAL)
  
  tumour_gt_df = current_vcf@gt
  colnames(tumour_gt_df)[2] = "TUMOR"
  
  tumour_info = getINFO(current_vcf) %>%
    as.data.frame() %>% 
    rename(INFO = ".")
  
  tumour_df = getFIX(current_vcf) %>%
    as.data.frame() %>% 
    mutate(ID = paste(samp, CHROM, POS, REF, sep = "_")) %>% 
    cbind(tumour_info, tumour_gt_df) %>% 
    left_join(normal_df) %>% 
    filter(!is.na(ALT))
  
  tumour_df_filt = tumour_df %>% 
    mutate(ID = paste(patient, CHROM, POS, REF, ALT, sep = "_")) %>% 
    filter(ID %in% gl_IDs)
  
  write_delim(tumour_df_filt, 
              file = file.path(maf_input_dir_blood,paste0(samp,"_vcf_filtered.txt")),
              delim = "\t",col_names = F)
}

