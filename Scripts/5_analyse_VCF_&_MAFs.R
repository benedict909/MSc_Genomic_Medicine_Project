# Script 5 - analyse MAF & VCF files

options(connectionObserver = NULL); library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
library(maftools) # BiocManager::install("PoisonAlien/maftools")
library(tidyverse)
library(grid)
library(venn)
library(ggtext) 


load2object <- function (filename){
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error! the file ", filename, 
            " was not found :("))
  NULL
}

# 0. define directories and load data ----
vcf_dir = "~/Masters_research_project/VCF_files"
per_sample_dir = "~/Masters_research_project/per_sample_vcf/"
biomarker_dir = "~/../Dropbox/Masters/9_Research_Project/Data/Biomarkers"
resdir = "~/../Dropbox/Masters/9_Research_Project/Report/Figures";  if(!dir.exists(resdir))dir.create(resdir)

load(file.path(vcf_dir,"vcf_LBsamples_filtered_allcolumns.RData")) # vcf_somatic
load(file.path(vcf_dir,"VCF_allsamples_unfiltered_allcolumns_annovar.RData"))
load(file.path(vcf_dir,"annot_table.RData")) # annot

plot_figs = c(TRUE, FALSE)[2]


# 1. Make annot table & define colours ----

# define colours for plotting 
all_cols = c("white","chartreuse3","orangered1","white",
             "white","#E454D6","goldenrod1", "white",
             "white",paste0(rep("turquoise",4),seq(1,4,1)))

names(all_cols) = c("**Responders**","Patient 1", "Patient 2", " ",
                    "**Non-responders**","Patient 3", "Patient 4","  ", 
                    "**Lung**",c("10ng","10ng1","20ng","40ng"))

pdac_cols = all_cols[!grepl("Lung|ng",names(all_cols))]

# n_weeks = c(0,6,12,26); names(n_weeks) =sort(unique(vcf_annotated$Timepoint))
# n_weeks_df = as.data.frame(n_weeks) %>%
#   rownames_to_column("Timepoint")
# 
# annot = vcf_annotated %>%
#   select(Sample, Patient, Timepoint, Tissue) %>%
#   distinct() %>%
#   mutate(Tissue = ifelse(Tissue == "U","Urine",Tissue),
#          Tissue = ifelse(Tissue == "PLASMA","Plasma",Tissue),
#          Treatment = ifelse(grepl("1|2",Patient), "Responder",NA),
#          Treatment = ifelse(grepl("3|4",Patient), "Non-responder",Treatment),
#          ) %>%
#   group_by(Sample) %>%
#   mutate(Amount = ifelse(grepl("14_P",Sample),unlist(strsplit(Sample,"_"))[3],NA),
#          graphcol = ifelse(!is.na(Amount),Amount,as.character(Patient)),
#          graphcol = gsub("_"," ",graphcol),
#          graphcol = factor(graphcol, levels = names(all_cols)),
#          pdaccol = ifelse(Patient!="Lung",as.character(graphcol),NA),
#          pdaccol = factor(pdaccol,levels = names(pdac_cols)),
#          Tumor_Sample_Barcode = Sample
#          ) %>%
#   left_join(n_weeks_df) %>%
#   arrange(Patient, Tissue, Timepoint)
# 
# save(annot,file =file.path(vcf_dir,"annot_table.RData"))


# 2. Sequencing stats ----

# mean average of Mb per sample passing sequencing quality control metrics in the gene panel
per_sample_dir = "~/Masters_research_project/per_sample_vcf/"; if(!dir.exists(per_sample_dir)) dir.create(per_sample_dir)
nbases = c()
mean_DP = c()
for(filename in list.files(per_sample_dir)){
  print(paste("Loading and adding", filename, "to the combined VCF"))
  
  current_vcf = load2object(file.path(per_sample_dir,filename))
  
  nbases_current = current_vcf %>% 
    filter(FILTER=="PASS") %>% 
    nrow()
  
  mean_DP_current = current_vcf %>% 
    filter(FILTER=="PASS") %>% 
    pull(DP) %>% 
    mean()

  rm(current_vcf) # to save memory
  names(nbases_current) = filename
  names(mean_DP_current) = filename
  
  if(filename == per_sample_files[1]){
    nbases = nbases_current
    mean_DP = mean_DP_current
  }else{
    nbases = c(nbases, nbases_current)
    mean_DP = c(mean_DP, mean_DP_current)
  }
}

# mean bases passing quality control (i.e. FILTER == "PASS")
mean(nbases) / 1000000
mean(nbases[!grepl("Blood",names(nbases))]) / 1000000
mean(nbases[grepl("Blood",names(nbases))]) / 1000000

# mean depth
mean(mean_DP) 
mean(mean_DP[!grepl("Blood",names(mean_DP))]) 
mean(mean_DP[grepl("Blood",names( mean_DP))]) 

# mean no. of blood variants
vcf_annotated %>% 
  filter(Tissue=="Blood"&Patient!="Lung", FILTER=="PASS") %>% 
  count(Sample) %>% 
  pull(n) %>% summary()

# mean no. of liquid biopsy variants
vcf_annotated %>% 
  filter(Tissue!="Blood"&Patient!="Lung", FILTER=="PASS") %>% 
  count(Sample) %>% 
  pull(n) %>% summary()

# ctDNA variant stats
vcf_somatic %>% 
  count(Sample) %>% 
  pull(n) %>% summary()

# restore variants stats
restore_counts = vcf_annotated %>% 
  filter(FILTER == "LowSupport" & Tissue != "Blood" & 
         Sample_ID %in% vcf_somatic$Sample_ID & alt_DP >= 3) %>%  # keep lowsupport variants with at least 3 ALT observations
  count(Sample,Patient, name = "restore_count") 

restore_counts %>% 
  pull(restore_count) %>% summary()

# corr
# germline variants
vcf_blood = vcf_annotated %>% 
  filter(Tissue == "Blood" & FILTER == "PASS") 

# PDAC somatic variants = variants not found in blood
blood_IDS = vcf_annotated %>% 
  filter(Tissue == "Blood" & FILTER == "PASS") %>% 
  pull(Sample_ID)

count_tab = vcf_annotated %>% 
  filter(FILTER == "PASS" & Tissue != "Blood" & !Sample_ID %in% blood_IDS, Patient!="Lung") %>% 
  count(Patient, name = "pass_count") %>% 
  left_join(restore_counts) 

count_tab %>% 
  ggplot(aes(x = pass_count, y = restore_count)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = T)
  
cor.test(count_tab$pass_count,count_tab$restore_count)

lm(count_tab$restore_count~count_tab$pass_count)

# somatic stats
table(vcf_somatic$Type)
nrow(vcf_somatic)

# gDNA/cfDNA comparison
table(blood_IDS %in% vcf_annotated$Sample_ID[vcf_annotated$Tissue!="Blood"])
length(unique(blood_IDS))


# 2. Tumour mutation burden (TMB) analysis ----

panel_df = data.frame("Sample" = unique(vcf_somatic$Sample), "n_megabases" = NA)

if(!file.exists(file.path(vcf_dir,"panel_df.R"))){
  for(mysamp in panel_df$Sample){
    current_vcf = load2object(file.path(per_sample_dir,paste0(mysamp,"_vcf.RData")))
    n_bases = current_vcf %>% 
      filter(FILTER == "PASS") %>% 
      nrow()
    panel_df$n_megabases[panel_df$Sample == mysamp] = n_bases / 1000000
  }
  save(panel_df, file = file.path(vcf_dir,"panel_df.R"))
}else{
  load(file.path(vcf_dir,"panel_df.R"))
}

# create input df
tmb_df = vcf_somatic %>% 
  count(Sample, name = "n_mut") %>% # count no. of mut per sample
  left_join(annot) %>% 
  left_join(panel_df) %>% 
  mutate(TMB = n_mut / n_megabases, # divide mut no. by panel size to calc TMB
         Patient_tissue = paste(Patient, Tissue, sep = "_")) %>% 
  group_by(graphcol, n_weeks) %>% 
  mutate(TMB_mean = mean(TMB),
         Type  = ifelse(Patient == "Lung","Lung","PDAC"),
         n_weeks = ifelse(Patient == "Lung",0,n_weeks))

# plot
if(plot_figs) png(file.path(resdir,"TMB_nolung.png"), width = 1500,height = 1000,res = 350)
ggplot(filter(tmb_df,Patient != "Lung")) + 
  geom_line( aes(x = n_weeks, y = TMB_mean, group = Patient)) + 
  geom_point(aes(x = n_weeks, y = TMB_mean, colour = pdaccol), size = 2.5) + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
  theme_bw() + theme(legend.title = element_blank(), legend.text = element_markdown()) + 
  scale_color_manual(drop=F, values = pdac_cols)+
  labs(y = "Tumour mutation burden (mut/Mb)", x = "Weeks after diagnosis")
if(plot_figs) dev.off()

p1_vals = tmb_df %>% 
  filter(Patient=="Patient_1",Timepoint!="T3") %>% 
  pull(TMB_mean)
(max(p1_vals) - min(p1_vals))/12 # no. of muts gained per week in patient 1

p2_vals = tmb_df %>% 
  filter(Patient=="Patient_2") %>% 
  pull(TMB_mean)
(max(p2_vals) - min(p2_vals))/26 # no. of muts gained per week in patient 1

# plot with lung
tmb_df$n_weeks[tmb_df$Patient=="Lung"] = 30 #so x-axis label is not diagnosis

tmb_plot = ggplot() +
  geom_point(data = filter(tmb_df, Patient == "Lung"),aes(x = n_weeks, y = TMB, colour = graphcol), size = 2.5) +
  geom_line(data =  filter(tmb_df, Patient != "Lung"), aes(x = n_weeks, y = TMB_mean, group = Patient)) + 
  geom_point(data = filter(tmb_df, Patient != "Lung"),aes(x = n_weeks, y = TMB_mean, colour = graphcol), size = 2.5) +
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,"")) + 
  theme_bw() + theme(legend.title = element_blank(), legend.text = element_markdown()) + 
  scale_color_manual(drop=F, values = all_cols)+
  labs(y = "Tumour mutation burden (mut/Mb)", x = "Weeks after diagnosis") +
  facet_grid(~ Type, scales = "free_x", space = "free")

tmb_gt = ggplot_gtable(ggplot_build(tmb_plot))
# gtable::gtable_show_layout(gt)
tmb_gt$widths[5] = 50*tmb_gt$widths[5] # increase the width of the lung facet

if(plot_figs) png(file.path(resdir,"TMB_lung.png"), width = 1500,height = 1000,res = 300)
grid.draw(tmb_gt)
if(plot_figs) dev.off()


# exonic TMB 
vcf_somatic %>% 
  filter(Func.refGene=="exonic") %>% 
  count(Sample, name = "n_mut") %>% 
  left_join(annot) %>% 
  group_by(Patient, n_weeks) %>% 
  mutate(TMB = n_mut / 1.2,
         TMB_mean = mean(TMB),
         Type  = ifelse(Patient == "Lung","Lung","PDAC"),
         n_weeks = ifelse(Patient == "Lung",0,n_weeks)) %>% 
  ggplot() + 
  geom_line( aes(x = n_weeks, y = TMB_mean, group = Patient)) + 
  geom_point(aes(x = n_weeks, y = TMB_mean, colour = Patient), size = 2.5) + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
  scale_y_continuous(limits = c(0,12.5), expand = c(0,0))+
  theme_bw() + theme(legend.title = element_blank()) + 
  labs(y = "Tumour mutation burden (mut/Mb)", x = "Weeks after diagnosis")

# however, as gene panel only 500 genes, neoantigens may be generated by mut in other genes of the genome
# therefore better do calc somatic mutation rate of the whole panel, not just coding mutations
# as this gives a more accurate readout of the likely relative no. of neoantigens present in the sample

# comaprison 

exonic_tmb = vcf_somatic %>% 
  filter(Func.refGene=="exonic") %>% 
  count(Sample, name = "n_mut") %>% 
  left_join(annot) %>% 
  group_by(graphcol, n_weeks) %>% 
  mutate(TMB = n_mut / 1.2,
         exonic_tmb = mean(TMB)) %>% 
  ungroup() %>% 
  dplyr::select(Sample, exonic_tmb)

tmb_df = tmb_df %>% 
  left_join(exonic_tmb)

tmb_corr = cor.test(tmb_df$TMB_mean, tmb_df$exonic_tmb)
tmb_text = paste("atop(italic(r)==",round(tmb_corr$estimate,2),",p.val==",
                 round(tmb_corr$p.value,11),")") # ?plotmath
corcol =  "royalblue4"

if(plot_figs) png(file.path(resdir,"TMB_correltaion.png"), width = 1500,height = 1000,res = 350)
# if(plot_figs) png(file.path(resdir,"TMB_correltaion_annot.png"), width = 1500,height = 1500,res = 350)
ggplot(tmb_df, aes(y = TMB_mean, x = exonic_tmb)) + 
  geom_point(aes(color = graphcol), size = 2.5) + 
  geom_smooth(method = "lm",se = F, colour = corcol, linetype = "dashed") + 
  scale_color_manual(drop=F, values = all_cols) + 
  theme_bw() + theme(legend.title = element_blank(), legend.text = element_markdown()) + 
  labs(y = "Total TMB (mut/Mb)", x = "Coding TMB (mut/Mb)") + 
  geom_label(aes(label = tmb_text, x = 8.5, y = 4), colour = corcol, fill = "white",
             size = 3, parse = T) + 
  scale_x_continuous(breaks = seq(0,14,2), limits = c(0,12.5), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(0,25,5),limits = c(0,27), expand = c(0,0))+ #limits = c(0,2), expand = c(0,0)
if(plot_figs) dev.off()

# corr with non-synomyomous SNV rate
nonsyn_tmb = vcf_somatic %>% 
  filter(ExonicFunc.refGene=="nonsynonymous SNV") %>% 
  count(Sample, name = "n_mut") %>% 
  left_join(annot) %>% 
  group_by(graphcol, n_weeks) %>% 
  mutate(TMB = n_mut / 1.2,
         nonsyn_tmb = mean(TMB)) %>% 
  ungroup() %>% 
  dplyr::select(Sample, nonsyn_tmb)

tmb_df = tmb_df %>% 
  left_join(nonsyn_tmb)

cor.test(tmb_df$TMB_mean, tmb_df$nonsyn_tmb)

# 3. Venn diagrams ----

IDs_list = c() # generate list of mutation IDs present in each sample 
for(i in 1:length(unique(vcf_somatic$Sample))){
  samp = unique(vcf_somatic$Sample)[i]
  current_IDs = vcf_somatic %>% 
    filter(Sample == samp) %>% 
    pull(ID)
  
  IDs_list = c(IDs_list, list(current_IDs))
  names(IDs_list)[i] = samp 
}

# patient 1
venn(x = IDs_list[annot$Sample[annot$Patient == "Patient_1" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F, ggplot = T) +
  ggtitle("Patient 1, all timepoints") + theme(plot.title = element_text(size = 20))

Patient_1_urine = unlist(IDs_list[annot$Sample[annot$Patient == "Patient_1" & annot$Tissue == "Urine"]])
Patient_1_plasma = unlist(IDs_list[annot$Sample[annot$Patient == "Patient_1" & annot$Tissue == "Plasma"]])

venn(x = list("Urine" = Patient_1_urine, "Plasma" = Patient_1_plasma),
     zcolor = "style", box = F, ggplot = T) +
  ggtitle("Patient 1, Urine vs plasma (all timepoints)") + theme(plot.title = element_text(size = 20))

# patient 2
venn(x = IDs_list[annot$Sample[annot$Patient == "Patient_2" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F,ggplot = T)+ 
  ggtitle("Patient 2") + theme(plot.title = element_text(size = 20))

#patient 3
venn(x = IDs_list[annot$Sample[annot$Patient == "Patient_3" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F,ggplot = T)+ 
  ggtitle("Patient 3") + theme(plot.title = element_text(size = 20))

# patient 4
venn(x = IDs_list[annot$Sample[annot$Patient == "Patient_4" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F,ggplot = T)+ 
  ggtitle("Patient 4") + theme(plot.title = element_text(size = 20))

# lung
venn(x = IDs_list[annot$Sample[annot$Patient == "Lung" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F,ggplot = T)+ 
  ggtitle("Lung") + theme(plot.title = element_text(size = 20))


# 4. MAF load & summary ----

# A. load & combine germline MAFs
maf_input_dir = "~/Masters_research_project/MAF_files"
input_files = list.files(maf_input_dir)[grepl("Blood", list.files(maf_input_dir))]

for(input_file in input_files){
  if(input_file == input_files[1]){
    maf_germline = read.maf(file.path(maf_input_dir, input_file),verbose = F, clinicalData = annot)
  }else{
    current_maf = read.maf(file.path(maf_input_dir, input_file), verbose = F, clinicalData = annot)
    maf_germline = merge_mafs(list(maf_germline, current_maf),verbose = F)
  }
}

maf_germline@data %>% 
  filter(FILTER == "PASS") %>% 
  left_join(annot) %>% 
  select(Patient,Hugo_Symbol,HGVSp, Consequence, CLIN_SIG, gnomAD_AF, SIFT, PolyPhen) %>% 
  arrange(Patient)

gl_keep = maf_germline@data %>% 
  filter(FILTER == "PASS", CLIN_SIG != "" | Hugo_Symbol == "MSH3") %>% # add MSH3 due to association with MSI
  left_join(annot) %>% 
  filter(Patient != "Lung") %>% 
  select(Patient,Hugo_Symbol,POS = Start_Position, REF = Reference_Allele,ALT = Tumor_Seq_Allele2, HGVSp, Consequence, CLIN_SIG, gnomAD_AF, SIFT, PolyPhen) %>% 
  arrange(Patient)

gl_keep

# write germline table 
gl_write = maf_germline@data %>% 
  filter(FILTER == "PASS", CLIN_SIG != "" | Hugo_Symbol == "MSH3") %>% # add MSH3 due to association with MSI
  left_join(annot) %>% 
  filter(Patient != "Lung") %>% 
  select(`gnomAD freq.` = gnomAD_AF, everything()) %>% 
  mutate(ID = paste0(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_")) %>% 
  left_join(vcf_somatic) %>% 
  select(Patient,Gene = Hugo_Symbol, `Transcript ID` = RefSeq, HGVSc, HGVSp, `Variant consequence` = Consequence, 
         CLIN_SIG, `gnomAD freq.`, `1000G`, SIFT, PolyPhen) %>% 
  arrange(Patient) %>% 
  mutate(Patient = gsub("Patient_","",Patient), 
         `Variant consequence` = gsub("_"," ",`Variant consequence`),
         PolyPhen = gsub("_"," ", PolyPhen),
         SIFT = gsub("_"," ", SIFT),
         CLIN_SIG = gsub("_"," ", CLIN_SIG))

write_delim(gl_write, 
            file = file.path(resdir,("germline_variants.txt")),
            delim = "\t",col_names = T)

# max gnomAD freq not designated as "common variant" by VEP 
maf_germline@data %>% 
  filter(!is.na(gnomAD_AF), FILTER != "common_variant") %>% 
  pull(gnomAD_AF) %>% 
  max()

# make df of germline muts to add to somatic MAFs as a "cnTable"
# gl_keep_genes = c("MLH1","MSH2","MSH3","B2M") # genes that are marked PASS in maf_germline@data$FILTER that are relevant to PDAC
gl_keep_genes = gl_keep$Hugo_Symbol 

gl_df = gl_keep %>% 
  left_join(annot) %>% 
  select(SYMBOL = Hugo_Symbol, Tumor_Sample_Barcode, Consequence) %>% 
  mutate(Consequence = ifelse(Consequence == "missense_variant", "GL_Missense", "GL_Startloss"))

adddf = as.data.frame(t(c("remove","CIRCB33_PLASMA","Variant"))); names(adddf) = names(gl_df)
gl_df = rbind(gl_df,adddf)

# B. Load somatic MAFs 

maf_input_dir = "~/Masters_research_project/MAF_files"
input_files = list.files(maf_input_dir)[!grepl("Blood|14", list.files(maf_input_dir))]

blacklist = c("CIRCB33_PLASMA.maf") # as only silent variants
for(input_file in input_files){
  if(!input_file %in% blacklist){  
    if(input_file == input_files[1]){
      maf_summary = read.maf(file.path(maf_input_dir, input_file), cnTable = c(), 
                              clinicalData = annot, verbose = F)
    }else{
      current_maf = read.maf(file.path(maf_input_dir, input_file), verbose = F, 
                             clinicalData = annot)
      maf_summary = merge_mafs(list(maf_summary, current_maf),verbose = F,cnTable = c())
    }
  }
}

# C. MAF summary 
# head(maf_summary@data, 2) # MAF similar to a VCF
getSampleSummary(maf_summary)
genesum = getGeneSummary(maf_summary); genesum #Shows gene summary.
getClinicalData(maf_summary) #shows clinical data associated with samples

variant_cols = c("#33A02C","#1f78b4","#000000","#535c68","#efa9ae","#f4e87c","white")
names(variant_cols) = c("Missense_Mutation","Frame_Shift_Del","Multi_Hit","pathway","GL_Missense",      
                        "GL_Startloss","Variant")

maf_summary = subsetMaf(maf_summary, query = "!grepl('Blood|14|_B', Tumor_Sample_Barcode)")

if(plot_figs) png(file.path(resdir,"variant_summary.png"), width = 2500,height = 1500,res = 300)
plotmafSummary(maf = maf_summary, rmOutlier = F, addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE, color = variant_cols, showBarcodes = F)
if(plot_figs) dev.off()

table(vcf_somatic$ExonicFunc.refGene[vcf_somatic$Patient!="Lung"])


# MAF colours 
patient_cols  = c("chartreuse3", "orangered1", "#E454D6","goldenrod1"); names(patient_cols) = paste0(rep("Patient_",4),c(1:4))
tissue_cols = c("tomato4", "goldenrod1"); names(tissue_cols) = c("Plasma", "Urine")
timepoint_cols = c("grey90","grey65", "grey40", "grey15"); names(timepoint_cols) = c( paste0(rep("T",4),c(0:3)))
response_cols = c("powderblue", "purple4"); names(response_cols) = c("Responder","Non-responder")
annot_cols = list(Patient = patient_cols, Timepoint = timepoint_cols, Tissue = tissue_cols, Treatment = response_cols)

# make annotation matrix for variant summary
wafcols = annot_cols
names(wafcols) = NULL

colourdf = as.data.frame(unlist(wafcols)) %>% 
  rownames_to_column("value") %>% 
  # mutate(value = gsub("Timepoint.|Patient.|Tissue.|Treatment.","",value),
  #        value = ifelse(nchar(value)==1,paste0("Patient_",value),value)) %>%
  select(value, colour = `unlist(wafcols)`)

waf_df = getSampleSummary(maf_summary) %>% 
  as.data.frame() %>% 
  select(Tumor_Sample_Barcode,total) %>% 
  left_join(annot) %>% 
  filter(Tissue!="Blood") %>% 
  select(Patient, Timepoint,Tissue,Treatment) %>% 
  mutate(x = 1:16) %>% 
  reshape2::melt(id.vars=c("x")) %>% 
  left_join(colourdf) %>% 
  mutate(y = rep(1:4,each=16))

if(plot_figs) png(file.path(resdir,"variant_summary_annot.png"), width = 2500,height = 1500,res = 300)
ggplot(waf_df, aes(x = x, y = y,fill = value)) + 
  geom_tile(color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
  scale_fill_manual(values = unlist(wafcols)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()


# 6. MAF Oncoplots ----

maf_input_dir = "~/Masters_research_project/MAF_files"
input_files = list.files(maf_input_dir)[!grepl("Blood|14", list.files(maf_input_dir))]

blacklist = c("CIRCB33_PLASMA.maf") # as only silent variants
for(input_file in input_files){
  if(!input_file %in% blacklist){  
    if(input_file == input_files[1]){
      maf_combined = read.maf(file.path(maf_input_dir, input_file), cnTable = gl_df, 
                              clinicalData = annot, verbose = F)
    }else{
      current_maf = read.maf(file.path(maf_input_dir, input_file), verbose = F, 
                             clinicalData = annot)
      maf_combined = merge_mafs(list(maf_combined, current_maf),verbose = F,cnTable = gl_df)
    }
  }
}
samp_order = annot %>% # sample order for oncoplot
  arrange(Patient, Tissue, Timepoint) %>% 
  pull(Tumor_Sample_Barcode)

# get KEGG pathways
somatic_genes = genesum %>% 
  filter(Frame_Shift_Del > 0 | Missense_Mutation > 0) %>% 
  pull(Hugo_Symbol)

oncoplot_genes = c(somatic_genes,gl_keep_genes)

# get human KEGG pathways
kegg = gage::kegg.gsets("hsa") # BiocManager::install("gage")

pathways_df = kegg$kg.sets %>% 
  reshape2::melt() %>% 
  dplyr::select(Entrez_ID = value, Pathway = L1) %>% 
  mutate(Entrez_ID = as.character(Entrez_ID))

# extract names of all genes in VCF in dataset
all_gene_names = vcf_annotated %>% 
  dplyr::select(Gene = Gene.refGene) %>% 
  separate_rows(Gene, sep = ";") %>% 
  distinct() %>% 
  pull(Gene)

# create a conversion matrix from entrez IDs to HUGO symbols
gene_convert = AnnotationDbi::mapIds(org.Hs.eg.db, all_gene_names, 'ENTREZID', 'SYMBOL') %>% 
  unlist() %>% 
  reshape2::melt() %>% as.data.frame() %>% 
  dplyr::select(Entrez_ID = "value") %>%
  rownames_to_column("Gene") %>% 
  left_join(pathways_df)

# make a dataframe of KEGG pathways
oncoplot_pathways = NA
for(mygene in oncoplot_genes){
  mydf = gene_convert %>% 
    filter(Gene == mygene) %>% 
    # dplyr::select(Gene, Entrez_ID) %>% 
    distinct()
  
  current_pathways = gene_convert %>% 
    filter(Gene == mygene)
  
  # add annots
  is_cancergene = sort(unique(ifelse(grepl("Cancer", current_pathways$Pathway, ignore.case = T),T,F)),
                       decreasing = T)[1]
  is_signaling_gene = sort(unique(ifelse(grepl("signaling pathway", current_pathways$Pathway),T,F)),
                           decreasing = T)[1]
  is_cellcycle_gene = sort(unique(ifelse(grepl("Cell cycle", current_pathways$Pathway),T,F)),
                           decreasing = T)[1]
  
  mydf = mydf %>% 
    mutate(Cancer_gene = ifelse(is_cancergene, T, F),
           signaling_gene = ifelse(is_signaling_gene, T, F),
           cell_cyle_gene = ifelse(is_cellcycle_gene, T, F))
  
  if(all(is.na(oncoplot_pathways))){
    oncoplot_pathways = mydf
  }else{
    oncoplot_pathways = rbind(oncoplot_pathways, mydf)
  }
}

# pathway dataframes
oncoplot_pathways %>% filter(grepl("pancreatic", Pathway, ignore.case = T))
oncoplot_pathways %>% filter(grepl("signaling", Pathway, ignore.case = T)) %>% head(4)
oncoplot_pathways %>% filter(grepl("cycle", Pathway, ignore.case = T))

# extract genes in pathways of interest from kegg df 
extract_genes = function(mypattern, pathways_df = oncoplot_pathways, filter_genes = NA){
  res_genes = pathways_df %>% 
    filter(grepl(mypattern, Pathway, ignore.case = T)) %>% 
    select(Gene) %>%
    distinct() %>%
    pull(Gene)
  res_genes = res_genes[!res_genes %in% filter_genes]
  return(res_genes)
}

cancer_genes = c(extract_genes("cancer|carcinoma"),"BRD4", "ZNF217", "FANCA", "MDC1", "SETBP1","PRKDC") # add genes manually from online searches
cycle_genes =  c(extract_genes("cycle",filter_genes = c(cancer_genes)), "CDK12")
signal_genes = c(extract_genes("signaling", filter_genes = c(cancer_genes)), "EGFL7",
                 "ABL2", "ZFHX3", "NF1", "INSR", "RNF43", "NOTCH1", "RECQL4", "PDGFRA")
immune_genes = c("PRDM1","B2M")

# combine into dataframe of pathways and pathway annotations for genes in oncoplot
pathways_2col = list("Cancer genes" = cancer_genes, "Cell cycle" = cycle_genes, 
                     "Cell signalling" = signal_genes, "Immune system"= immune_genes) %>% 
  reshape2::melt() %>% 
  rename(Gene = value, Pathway = L1)

# make right hand side bar
gl_muts = gl_df %>% left_join(annot) %>% select(Gene = SYMBOL, Patient) %>% distinct()

# count of no. of patients with mutations in each pathway
pathways_count = maf_combined@data %>% 
  select(where(~!all(is.na(.x)))) %>% 
  left_join(annot) %>%
  select(Gene = SYMBOL, Patient) %>% 
  rbind(gl_muts) %>% 
  left_join(pathways_2col) %>% 
  select(Patient, SYMBOL = Pathway) %>% 
  distinct() %>% 
  count(SYMBOL, name = "No. of patients") 

# count of no. of germline muts per gene per patient
gl_count = gl_df %>% 
  left_join(annot) %>%
  select(SYMBOL, Patient) %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  count(name = "No. of patients")

# combine pathway & germline counts with somatic counts
oncoplot_rightbar = maf_combined@data %>% 
  select(where(~!all(is.na(.x)))) %>% 
  left_join(annot) %>%
  select(SYMBOL, Patient) %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  count(name = "No. of patients") %>%
  as.data.frame() %>% 
  rbind(pathways_count, gl_count) %>% 
  filter(!is.na(SYMBOL),
         !SYMBOL %in% SYMBOL[duplicated(.$SYMBOL)]) 
# plot 
clin_feats = c("Patient","Timepoint" ,"Tissue","Treatment") # annots to include under oncoplot
names(annot_cols$Patient) = gsub(" ","_",names(annot_cols$Patient))

maf_onco = subsetMaf(maf_combined, query = "!grepl('Blood|14', Tumor_Sample_Barcode)")


if(plot_figs) png(file.path(resdir,"Oncoplot_PDAC.png"), width = 3000,height = 3000,res = 300)
oncoplot(maf = maf_onco, showTumorSampleBarcodes = F,
         cohortSize = 17,clinicalFeatures = clin_feats, sampleOrder = samp_order, 
         annotationColor = annot_cols, colors = variant_cols, writeMatrix = F, drawBox = F, 
         fontSize = 0.75, anno_height = 1.2, gene_mar = 8, drawRowBar = T, annotationFontSize = 1,
         legendFontSize = 1, showTitle = F, #genesToIgnore = "remove",
         pathways = pathways_2col, rightBarData = oncoplot_rightbar, rightBarLims = c(0,4)
)
if(plot_figs) dev.off()

pathways_2col[nrow(pathways_2col)+1,] = c("remove","other") # add variant text to graph 
if(plot_figs) png(file.path(resdir,"Oncoplot_PDAC_annot.png"), width = 3000,height = 3000,res = 300)
oncoplot(maf = maf_onco, showTumorSampleBarcodes = F,
         cohortSize = 17,clinicalFeatures = clin_feats, sampleOrder = samp_order, 
         annotationColor = annot_cols, colors = variant_cols, writeMatrix = F, drawBox = F, 
         fontSize = 0.75, anno_height = 1.2, gene_mar = 8, drawRowBar = T, annotationFontSize = 1,
         legendFontSize = 1, showTitle = F, genesToIgnore = "remove",
         pathways = pathways_2col, rightBarData = oncoplot_rightbar, rightBarLims = c(0,4)
)
if(plot_figs) dev.off()



# alt consistency - better version
mafannot = as.data.frame(maf_onco@data) %>% 
  dplyr::select(where(~!all(is.na(.x)))) %>% 
  filter(!grepl("GL",Variant_Classification)) %>% 
  left_join(annot) %>% 
  mutate(Patient_Tissue =  paste(Patient, Tissue, sep = "_"))

mutpersamp_df = mafannot %>% 
  group_by(Patient_Tissue, Hugo_Symbol) %>% 
  count(Hugo_Symbol, name = "muts_per_samp")

# no. of samples per patient
prevalence_df = annot %>% 
  filter(Tissue != "Blood", Patient != "Lung") %>% 
  mutate(Patient_Tissue =  paste(Patient, Tissue, sep = "_")) %>% 
  group_by(Patient_Tissue) %>% 
  count(Patient, name = "samps_per_patient") %>% 
  left_join(mutpersamp_df) %>% 
  mutate(all_samples = ifelse(samps_per_patient == muts_per_samp, "yes","no"))

table(prevalence_df$all_samples) # no. of sample/gene combos = sum of yes/no
# there are 48 somatic gene-sample combos in the oncoplot

# % of mutations that 1st appear in T0
mafannot %>% 
  filter(Timepoint == "T0") %>% 
  select(Patient_Tissue,Hugo_Symbol) %>% 
  nrow() / nrow(prevalence_df)

if(plot_figs) png(file.path(resdir,"Oncoplot_consistency.png"), width = 1000,height = 700,res = 300)
prevalence_df %>% 
  mutate(prevalence = ifelse(all_samples == "yes","All timepoints", "Some timepoints")) %>% 
  ggplot(aes(x = prevalence, fill = all_samples))+
    geom_bar(stat = "count", width = 0.5) + 
    scale_y_continuous(expand = c(0,0), limits  = c(0,27)) + 
    theme_bw()+
    labs(x = "Mutation presence in each sample", y = "No. of genes")+
    theme(legend.position = "none")
if(plot_figs) dev.off()


# 7. Mutant allele frequency (MAF) analysis ----

summary(vcf_somatic$gt_VF)

maf_df = vcf_somatic %>% 
  group_by(Sample) %>% 
  summarise(median_MAF = median(gt_VF)) %>% # calculate median MAF in each sample
  left_join(annot) %>% 
  group_by(Patient, Timepoint) %>% 
  mutate(mean_median_MAF = mean(median_MAF), # calc mean MAF as multiple samps/timepoint in Patient 1 
         Type  = ifelse(Patient == "Lung","Lung","PDAC"),
         n_weeks = ifelse(Patient == "Lung",0,n_weeks)) 

# Plot MAF
if(plot_figs) png(file.path(resdir,"MAF_nolung.png"), width = 1500,height = 1000,res = 350)
MAF_nolung = ggplot(filter(maf_df,Patient!="Lung"), aes(x = n_weeks, y = mean_median_MAF)) + 
  geom_line( aes(group = Patient)) + 
  geom_point(aes(colour = pdaccol), size = 2.5) + 
  scale_color_manual(drop=F, values = pdac_cols)+
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
  scale_y_continuous(limits = c(0,.5), expand = c(0,0))+
  labs(y = "Median mutation allele frequency", x = "Weeks after diagnosis") + 
  theme_bw() + theme(legend.title = element_blank(),legend.text =  element_markdown(),
                     legend.position = "bottom")
if(plot_figs) dev.off()


# Plot MAF with lung
maf_df$n_weeks[maf_df$Patient=="Lung"] = 30 #so x-axis label is not diagnosis

maf_plot = ggplot() +
  geom_jitter(data =  filter(maf_df, Patient == "Lung"),aes(colour = graphcol, x = n_weeks, y = median_MAF),
              size = 2.5,height = 0.01, width = 0) + 
  geom_line(data =   filter(maf_df, Patient != "Lung"), aes(group = Patient, x = n_weeks, y = mean_median_MAF)) +
  geom_point(data =  filter(maf_df, Patient != "Lung"),aes(colour = graphcol, x = n_weeks, y = mean_median_MAF), size = 2.5) +
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,"")) +
  scale_y_continuous(limits = c(0,.5), expand = c(0,0))+
  theme_bw() + theme(legend.title = element_blank(), legend.text =  element_markdown()) + 
  scale_color_manual(drop=F, values = all_cols)+
  labs(y = "Median mutation allele frequency", x = "Weeks after diagnosis") + 
  facet_grid(~ Type, scales = "free_x", space = "free")

maf_gt = ggplot_gtable(ggplot_build(maf_plot))
# gtable::gtable_show_layout(gt)
maf_gt$widths[5] = 50*maf_gt$widths[5] # increase the width of the lung facet

if(plot_figs) png(file.path(resdir,"MAF_lung.png"), width = 1500,height = 1000,res = 300)
grid.draw(maf_gt)
if(plot_figs) dev.off()


# calculate % changes
# a. % decrease in MAF in Patient 1 between T0 and T1
patient_1_T0_MAF = unique(maf_df$mean_median_MAF[maf_df$Patient == "Patient_1" & maf_df$Timepoint == "T0"])
patient_1_T1_MAF = unique(maf_df$mean_median_MAF[maf_df$Patient == "Patient_1" & maf_df$Timepoint == "T1"])
(patient_1_T0_MAF - patient_1_T1_MAF) * 100 / patient_1_T0_MAF 

patient_1_T2_MAF = unique(maf_df$mean_median_MAF[maf_df$Patient == "Patient_1" & maf_df$Timepoint == "T2"])
(patient_1_T0_MAF - patient_1_T2_MAF) * 100 / patient_1_T0_MAF 

# b. % increase in MAF in Patient 4 between T0 and T1
patient_4_T0_MAF = unique(maf_df$mean_median_MAF[maf_df$Patient == "Patient_4" & maf_df$Timepoint == "T0"])
patient_4_T1_MAF = unique(maf_df$mean_median_MAF[maf_df$Patient == "Patient_4" & maf_df$Timepoint == "T1"])
(patient_4_T1_MAF - patient_4_T0_MAF) * 100 / patient_4_T0_MAF

# plot patient 1 pathways MAF 

# maf_combo_df = as.data.frame(maf_combined@data) %>% 
#   dplyr::select(where(~!all(is.na(.x)))) %>% 
#   filter(Hugo_Symbol %in% somatic_genes) %>% 
#   dplyr::select(Sample = Tumor_Sample_Barcode, CHROM = Chromosome, POS = Start_Position, 
#                 REF = Reference_Allele, ALT = Tumor_Seq_Allele2) %>% 
#   left_join(annot) %>%
#   mutate(Sample_ID = paste(Sample, CHROM, POS, REF, ALT, sep = "_"),
#          Patient_ID = paste(Patient, CHROM, POS, REF, ALT, sep = "_"))
#          
# driver_IDs = unique(maf_combo_df$Patient_ID)

pathway_cols = c("white","skyblue2","orchid4","seagreen2","yellow3")
names(pathway_cols) = c("**Pathway**",unique(pathways_2col$Pathway[pathways_2col$Pathway!="other"]))

maf_p1_pathways_df = vcf_somatic %>% 
  mutate(Tissue = ifelse(Tissue == "U", "Urine","Plasma")) %>% 
  filter(Gene.refGene %in% oncoplot_genes, Func.refGene=="exonic") %>% 
  left_join(annot) %>% 
  rename(Gene = Gene.refGene) %>% 
  group_by(Patient, Gene, Timepoint) %>% 
  mutate(mean_VF = median(gt_VF),
         mean_DP = mean(alt_DP)) %>% 
  ungroup() %>% 
  left_join(pathways_2col) %>% 
  filter(Patient == "Patient_1", !is.na(Pathway)) %>% 
  group_by(Timepoint, Pathway) %>% 
  mutate(mean_VF_allgenes = mean(gt_VF),
         Pathway = factor(Pathway, levels = names(pathway_cols))) 

maf_p1_pathways_plot = maf_p1_pathways_df %>% 
  ggplot(aes(x = n_weeks, y = mean_VF_allgenes)) +
    geom_line( aes(group = Pathway), colour = "black") + 
    geom_point(size=2.5, aes(colour=Pathway)) + theme_bw() + 
    labs(x = "Weeks after diagnosis",y = "Median mutation allele frequency") + 
    scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
    scale_y_continuous(limits = c(0,.5), expand = c(0,0))+
    scale_colour_manual(values = pathway_cols, drop = F) + 
    # guides(colour = guide_legend(nrow = 2)) + 
    theme(legend.text =element_markdown(), legend.title = element_blank()) 

if(plot_figs) png(file.path(resdir,"MAF_patient1_pathways.png"), width = 1500,height = 1000,res = 350)
maf_p1_pathways_plot
if(plot_figs) dev.off()

# if(plot_figs) png(file.path(resdir,"MAF_bothpanels.png"), width = 3000,height = 1200,res = 350)
# ggpubr::ggarrange(MAF_nolung, maf_p1_pathways)
# if(plot_figs) dev.off()

maf_p1_pathways_df %>% 
  select(Pathway, Timepoint, mean_VF_allgenes) %>% 
  distinct() %>% 
  filter(grepl("Cancer",Pathway))

patient_1_T0_MAF = unique(maf_p1_pathways_df$mean_VF_allgenes[grepl("Cancer",maf_p1_pathways_df$Pathway) & maf_p1_pathways_df$Timepoint == "T0"])
patient_1_T1_MAF = unique(maf_p1_pathways_df$mean_VF_allgenes[grepl("Cancer",maf_p1_pathways_df$Pathway) & maf_p1_pathways_df$Timepoint == "T1"])
patient_1_T2_MAF = unique(maf_p1_pathways_df$mean_VF_allgenes[grepl("Cancer",maf_p1_pathways_df$Pathway) & maf_p1_pathways_df$Timepoint == "T2"])

(patient_1_T0_MAF - patient_1_T1_MAF) * 100 / patient_1_T0_MAF 
(patient_1_T0_MAF - patient_1_T2_MAF) * 100 / patient_1_T0_MAF 


# 7. Study driver genes from MAF ----
print_mut = function(gene_name = NA, input_df = as.data.frame(maf_summary@data), 
                     annot_df = annot) {
  if (all(is.na(gene_name))) 
    gene_name = unique(input_df$SYMBOL)
  returndf = input_df %>%
    dplyr::select(where(~!all(is.na(.x)))) %>%
    dplyr::select(Tumor_Sample_Barcode, SYMBOL, HGVSc,RefSeq, HGVSp, CLIN_SIG) %>%
    filter(SYMBOL %in% gene_name) %>%
    distinct() %>%
    left_join(annot_df) %>%
    arrange(Patient, SYMBOL, Timepoint)
  return(returndf)
}

print_mut(c("KRAS"))# G12* in KRAS most commonly mutated locus in PDAC wt KRAS = better prognosis

maf_sprint_mut(c("TP53"))
print_mut("PIK3R2")

# MLH1, MSH2 & MSH3 are DNA mismatch repair genes. Mutations in these genes are associated with microsatelite instability (MSI), a marker predictive if response to immunotherapy.
# ZFHX3 mutations are an independent predictive biomarker for NSCLC patients receiving ICI treatment (Zhang et al., 2021, Cancer Immunol Immunother.)
# Truncating variants in B2M (ß2-microglobulin), a component of MHC class I molecules, are associated with resistance to immunotherapy in melanoma (Zaretsky et al., 2016, NEJM).

maf_combined@data %>% 
  select(where(~!all(is.na(.x)))) %>% 
  filter(SYMBOL == "PRKDC") %>%  
  select(SYMBOL, HGVSc, HGVSp, PolyPhen, SIFT) %>% 
  distinct() # but not looking for pathogenicity, just change in structure that would mean response

p1_samps = vcf_somatic %>% filter(grepl("1",Patient)) %>% pull(Sample) %>% unique()
maf_p1 = subsetMaf(maf_combined, query = "Tumor_Sample_Barcode %in% p1_samps")

if(plot_figs) png(file.path(resdir,"PRKDC_lollipop_text.png"), width = 15000,height = 5000,res = 1000)
lollipopPlot(maf = maf_p1, gene = 'PRKDC',AACol = 'HGVSp', showMutationRate = TRUE, labelPos = "all",
             showDomainLabel = F, labelOnlyUniqueDoamins = T,)
if(plot_figs) dev.off()

if(plot_figs) png(file.path(resdir,"PRKDC_lollipop_image.png"), width = 1500,height = 700,res = 200)
lollipopPlot(maf = maf_p1, gene = 'PRKDC',AACol = 'HGVSp', showMutationRate = F,pointSize = 3,
             showDomainLabel = F, labelOnlyUniqueDoamins = T,titleSize = c(0,0),printCount = F)
if(plot_figs) dev.off()

rev(unique(maf_summary@data$HGVSc[maf_summary@data$Hugo_Symbol=="PRKDC"]))
unique(maf_summary@data$RefSeq[maf_summary@data$Hugo_Symbol=="PRKDC"])
rev(unique(maf_summary@data$HGVSp[maf_summary@data$Hugo_Symbol=="PRKDC"]))


# See Tan et al. (2020, Journal for ImmunoTherapy of Cancer) for more info on PRKDC mut association with immunotherapy


# 8. Biomarker analysis ----
biomarkerfiles = list.files(biomarker_dir)[!grepl(".zip",list.files(biomarker_dir))]

# load & combine biomarker .txt files
for(filename in biomarkerfiles){
  current_file = read.delim2(file.path(biomarker_dir,filename))
  current_file = current_file[6:nrow(current_file),]
  names(current_file) = c("Variable","Observation")
  current_file = current_file %>% 
    filter(Variable != "",Observation != "") %>% 
    mutate(Variable = gsub(":","",Variable),
           Variable = gsub(" ","_",Variable)) %>% 
    column_to_rownames("Variable") %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(Amount = NA)
  
  rownames(current_file) = NULL
  
  # convert character numbers to numeric
  numeric_cols = names(current_file)[6:13]
  current_file = current_file %>% 
    mutate(across(all_of(numeric_cols), as.numeric))
  
  # fix sample names 
  filename_split = unlist(strsplit(filename, "_"))
  if(length(filename_split) == 2)  filename_split = unlist(strsplit(filename_split, split = "(?<=[0-9])",
                                                                    perl = T)) # fix for when there is no underscore between sample and tissue
  if(!is.na(suppressWarnings(as.numeric(filename_split[2])))){ # fix for when there are 2 numbers in the sample name 
    filename_split[1] = paste0(filename_split[1],filename_split[2])
    filename_split[2] = filename_split[3]
  }
  samp_name = paste(filename_split[1],filename_split[2], sep = "_")
  lastchar = as.numeric(substr(samp_name,nchar(samp_name),nchar(samp_name)))
  if(!is.na(lastchar)) samp_name = substr(samp_name,1,nchar(samp_name)-1)
  
  curr_amount = NA
  if(grepl("ng", filename)) curr_amount = unlist(strsplit(unlist(strsplit(filename, "PLASMA"))[2],"_"))[1]
  
  current_file = current_file %>% 
    mutate(Amount = curr_amount,
           Sample_ID = samp_name,
           Filename = filename) %>% 
    dplyr::select(Sample = Sample_ID, Amount, Filename, Total_TMB, Nonsynonymous_TMB, 
                  Coding_Region_MB = Coding_Region_Size_in_Megabases, 
                  N_passing_eligble_vars = Number_of_Passing_Eligible_Variants,
                  N_passing_eligble_nonsyn_vars = Number_of_Passing_Eligible_Nonsynonymous_Variants, 
                  Usable_MSI_Sites,Total_Microsatellite_Sites_Unstable, Percent_Unstable_Sites)
  
  if(filename == biomarkerfiles[1]){
    biomarkers_df = current_file
  }else{
    biomarkers_df = rbind(biomarkers_df, current_file)
  }
}

biomarkers_df = biomarkers_df %>%
  mutate(Sample = ifelse(!is.na(Amount),paste(Sample,Amount,sep = "_"),Sample),
         Sample = ifelse(Sample == "CIRCB35_PLASMA","CIRCB35a_PLASMA",Sample)) %>% 
  left_join(annot) %>% 
  mutate(Type = ifelse(grepl("14",Sample),"Lung","PDAC"),
         n_weeks = ifelse(is.na(n_weeks),33,n_weeks)) %>% 
  group_by(Patient,Timepoint) %>% 
  mutate(Mean_perc_unstable = mean(Percent_Unstable_Sites))  

biomarkers_df %>% arrange(Patient, Tissue, Timepoint) %>% 
  dplyr::select(-c(N_passing_eligble_vars,N_passing_eligble_nonsyn_vars,
                   Tumor_Sample_Barcode, Filename))

# Plot Microsatellite instability (MSI)
# define function
plot_MSI = function(input, grouping = "Tissue", yvar= "Percent_Unstable_Sites", 
                    ylab = "Unstable microsatellite sites (%)",linecol="black"){
  if(length(unique(input$Patient))==1){
    mytitle = unique(input$Patient)
  }else{
    mytitle = "All patients"
  }
  if(linecol == "black"){
    ggplot(input, aes(x = n_weeks, y = !!as.name(yvar))) + 
      geom_line( aes(group = !!as.name(grouping)),colour=linecol) + 
      geom_point(aes(colour = !!as.name(grouping)), size = 2.5) + 
      scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
      labs(y = ylab, x = "Weeks after diagnosis", title = mytitle) + 
      theme_bw() + theme(legend.title = element_blank())
  }else{
    ggplot(input, aes(x = n_weeks, y = !!as.name(yvar))) + 
      geom_line( aes(group = !!as.name(grouping),colour=!!as.name(grouping))) + 
      geom_point(aes(colour = !!as.name(grouping)), size = 2.5) + 
      scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
      labs(y = ylab, x = "Weeks after diagnosis", title = mytitle) + 
      theme_bw() + theme(legend.title = element_blank())
  }
}

# patient 1
biomarkers_df %>% 
  filter(Patient == "Patient_1") %>% 
  plot_MSI()

biomarkers_df %>% 
  filter(Patient == "Patient_1", Tissue != "Blood") %>% 
  group_by(Timepoint) %>% 
  mutate(Mean_perc_unstable = mean(Percent_Unstable_Sites)) %>% 
  plot_MSI(yvar = "Mean_perc_unstable")

# Patient 2
biomarkers_df %>% 
  filter(Patient == "Patient_2") %>% 
  plot_MSI()

# Patient 3 
biomarkers_df %>% 
  filter(Patient == "Patient_3") %>% 
  plot_MSI()

# Patient 4
biomarkers_df %>% 
  filter(Patient == "Patient_4") %>% 
  plot_MSI()

# all patients
biomarkers_df %>% 
  filter(Tissue != "Blood") %>%
  group_by(Timepoint, Patient) %>% 
  mutate(Mean_perc_unstable = mean(Percent_Unstable_Sites)) %>%  
  plot_MSI(grouping = "Patient", yvar = "Mean_perc_unstable")

# calc corr of MSI increase over time
lmdf = biomarkers_df %>% 
  filter(Treatment=="Responder") %>% 
  select(n_weeks, Mean_perc_unstable) %>% 
  distinct()

responder_corr = cor.test(x=lmdf$Mean_perc_unstable,y=lmdf$n_weeks, method = "pearson")
cor_text = paste("atop(italic(r)==",round(responder_corr$estimate,2),",p.val==",
                 round(responder_corr$p.value,3),")") # ?plotmath
corcol =  "royalblue4"

cor_annot = data.frame(Type = "PDAC", x = 22, y = 1)

msiplot = ggplot() + 
  geom_line(data=filter(biomarkers_df,Patient != "Lung"),
            aes(group = Patient,x = n_weeks, y = Mean_perc_unstable)) +
  geom_point(data=filter(biomarkers_df,Patient == "Lung"),
             aes(colour = graphcol,x = n_weeks, y = Percent_Unstable_Sites), size = 2.5) + 
  geom_point(data=filter(biomarkers_df,Patient != "Lung"),
               aes(colour = graphcol,x = n_weeks, y = Mean_perc_unstable), size = 2.5) +
  scale_colour_manual(values = all_cols, drop=F) +
  theme_bw() + theme(legend.title = element_blank(), legend.text = element_markdown()) + 
  facet_grid(~ Type, scales = "free_x", space = "free") +
  scale_x_continuous(breaks = c(0,6,12,18,24,30,36), labels = c("Diagnosis",6,12,18,24,30,36)) +
  scale_y_continuous(limits = c(0,8.7), expand = c(0,0)) +
  labs(y = "Unstable microsatellite sites (%)", x = "Weeks after diagnosis")  +
  geom_smooth(data=filter(biomarkers_df,Treatment=="Responder"),method = "lm", se = FALSE, 
              colour = corcol, linetype = "dashed", size = 0.7, #show.legend=TRUE,
              aes(x = n_weeks, y = Mean_perc_unstable)) + 
  geom_label(data = cor_annot, aes(label = cor_text, x = 22, y = 1, group = Type), colour = corcol, fill = "white",
             size = 3, parse = T) 
  
msi_gt = ggplot_gtable(ggplot_build(msiplot))
# gtable::gtable_show_layout(msi_gt)
msi_gt$widths[5] = 50*msi_gt$widths[5] # increase the width of the lung facet

if(plot_figs) png(file.path(resdir,"MSI_lung.png"),  width = 1600,height = 1000,res = 300)
grid.draw(msi_gt)
if(plot_figs) dev.off()



if(plot_figs) png(file.path(resdir,"MSI_nolung.png"),  width = 1500,height = 1000,res = 300)
ggplot() + 
  geom_smooth(data=filter(biomarkers_df,Treatment=="Responder"),method = "lm", se = FALSE, 
              colour = corcol, linetype = "dashed", size = 0.7, #show.legend=TRUE,
              aes(x = n_weeks, y = Mean_perc_unstable)) +
  geom_line(data=filter(biomarkers_df,Patient != "Lung"),
            aes(group = Patient,x = n_weeks, y = Mean_perc_unstable)) +
  geom_point(data=filter(biomarkers_df,Patient != "Lung"), 
             aes(colour = pdaccol,x = n_weeks, y = Mean_perc_unstable), size = 2.5) +
  geom_label(aes(label = cor_text, x = 22, y = 1), colour = corcol, fill = "white",
             size = 3, parse = T) + 
  scale_colour_manual(values = pdac_cols, drop=F) +
  theme_bw() + theme(legend.title = element_blank(), legend.text = element_markdown()) + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) +
  scale_y_continuous(limits = c(0,7), expand = c(0,0)) +
  labs(y = "Unstable microsatellite sites (%)", x = "Weeks after diagnosis") 
if(plot_figs) dev.off()

# slope of lm
reg_output = summary(lm(formula= Mean_perc_unstable~n_weeks, data = lmdf))
reg_output$coefficients["n_weeks","Estimate"] # % MSI increase per week
reg_output$coefficients["(Intercept)","Estimate"]

sd(biomarkers_df$Percent_Unstable_Sites[biomarkers_df$Patient=="Lung"])

# Biomarkers file TMB
biomarkers_df = biomarkers_df %>% 
  filter(Tissue != "Blood") %>% 
  group_by(Timepoint, Patient) %>% 
  mutate(Mean_TMB = mean(Total_TMB))

tmbplot = ggplot() + 
  geom_line(data=filter(biomarkers_df,Patient != "Lung"), 
            aes(group = Patient,x = n_weeks, y = Mean_TMB)) +
  geom_point(data=filter(biomarkers_df,Patient == "Lung"),
             aes(colour = Patient,x = n_weeks, y = Total_TMB), size = 2.5) + 
  geom_point(data=filter(biomarkers_df,Patient != "Lung"),
             aes(colour = Patient,x = n_weeks, y = Mean_TMB), size = 2.5) + 
  theme_bw() + theme(legend.title = element_blank()) + 
  facet_grid(~ Type, scales = "free_x", space = "free") + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c(0,6,12,18,24,30)) +
  scale_y_continuous(limits = c(0,4.2), expand = c(0,0)) +
  labs(y = "TMB", x = "Weeks after diagnosis")  

tmb_gt2 = ggplot_gtable(ggplot_build(tmbplot))
# gtable::gtable_show_layout(gt)
tmb_gt2$widths[5] = 50*tmb_gt2$widths[5] # increase the width of the lung facet
grid.draw(tmb_gt2)

biomarkers_df %>% 
  filter(Tissue != "Blood") %>% 
  group_by(Timepoint, Patient) %>% 
  mutate(Mean_TMB = mean(Nonsynonymous_TMB)) %>%
  ggplot(aes(x = n_weeks, y = Mean_TMB)) + 
  geom_line( aes(group = Patient)) + 
  geom_point(aes(colour = Patient), size = 2.5) + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
  labs(y = "Nonsynonymous_TMB", x = "Weeks after diagnosis") + 
  theme_bw() + theme(legend.title = element_blank())

biomarkers_df %>% 
  filter(Tissue != "Blood") %>% 
  group_by(Timepoint, Patient) %>% 
  mutate(Mean_TMB = mean(Coding_Region_MB)) %>%
  ggplot(aes(x = n_weeks, y = Mean_TMB)) + 
  geom_line( aes(group = Patient)) + 
  geom_point(aes(colour = Patient), size = 2.5) + 
  scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
  labs(y = "Coding_Region_MB", x = "Weeks after diagnosis") + 
  theme_bw() + theme(legend.title = element_blank())

# 9. Lung analysis ----

vcf_lung = filter(vcf_somatic,Patient == "Lung")
vcf_lung_gl = vcf_annotated %>% 
  filter(Sample == "CIRCB14_Blood", !is.na(ALT), FILTER=="PASS")

vcf_lung %>% 
  count(Sample)

# all 8 somatic mutations in 10ng samples are shared
all(sort(vcf_lung$ID[vcf_lung$Sample=="CIRCB14_PLASMA_10ng"]) == sort(vcf_lung$ID[vcf_lung$Sample=="CIRCB14_PLASMA_10ng1"]))

# unique 3 in 20ng 
ID40ng = vcf_lung$ID[grepl("40ng",vcf_lung$Sample)]
ID20ng = vcf_lung$ID[grepl("20ng",vcf_lung$Sample)]

# Venn diagram for lung
IDs_list = c() # generate list of mutation IDs present in each sample 

for(i in 1:length(unique(vcf_lung$Sample))){
  samp = gsub("CIRCB14_PLASMA_","",unique(vcf_lung$Sample)[i])
  current_IDs = vcf_lung %>% 
    filter(gsub("CIRCB14_PLASMA_","",Sample) == samp) %>% 
    pull(ID)
  IDs_list = c(IDs_list, list(current_IDs))
  names(IDs_list)[i] = samp 
}
annot = annot %>% 
  mutate(Amount = ifelse(grepl("CIRCB14_PLASMA",Sample),gsub("CIRCB14_PLASMA_","",Sample),NA))

if(plot_figs) png(file.path(resdir,"lung_venn.png"),  width = 1000,height = 1000,res = 300)
venn(x = IDs_list[annot$Amount[annot$Patient == "Lung" & annot$Tissue != "Blood"]],
     zcolor = "style", box = F,ggplot = F, sncs = 1, ilcs = 1,ellipse = T)#+ 
  # ggtitle("Lung") + theme(plot.title = element_text(size = 20))
if(plot_figs) dev.off()


# oncoplots lung
# A. load & combine germline MAFs
maf_input_dir = "~/Masters_research_project/MAF_files"
input_files = list.files(maf_input_dir)[grepl("14_Blood", list.files(maf_input_dir))]

for(input_file in input_files){
  if(input_file == input_files[1]){
    maf_germline = read.maf(file.path(maf_input_dir, input_file),verbose = F, clinicalData = annot)
  }else{
    current_maf = read.maf(file.path(maf_input_dir, input_file), verbose = F, clinicalData = annot)
    maf_germline = merge_mafs(list(maf_germline, current_maf),verbose = F)
  }
}

gl_lung = maf_germline@data %>% 
  as.data.frame() %>% 
  filter(FILTER == "PASS", CLIN_SIG != "") %>% # add MSH3 due to association with MSI
  left_join(annot) %>% 
  filter(Patient == "Lung") %>% 
  select(Patient,Hugo_Symbol,HGVSp, Consequence, CLIN_SIG, gnomAD_AF, SIFT, PolyPhen) %>% 
  arrange(Patient)

gl_lung
# Rad50, a protein involved in DNA double-strand break repair.
# SUFU has also been found to have a crucial role in tumour suppression
gl_keep_lung = gl_lung$Hugo_Symbol 

gl_df_lung = gl_lung %>% 
  left_join(annot) %>% 
  select(SYMBOL = Hugo_Symbol, Tumor_Sample_Barcode, Consequence) %>% 
  mutate(Consequence = ifelse(Consequence == "missense_variant", "GL missense", "GL splice acceptor"))


# B. Load somatic MAFs 

maf_input_dir = "~/Masters_research_project/MAF_files"
input_files = list.files(maf_input_dir)[grepl("14_PLASMA", list.files(maf_input_dir))]

for(input_file in input_files){
  if(input_file == input_files[1]){
    maf_lung = read.maf(file.path(maf_input_dir, input_file), cnTable = gl_df_lung, 
                        clinicalData = annot, verbose = F)
  }else{
    current_maf = read.maf(file.path(maf_input_dir, input_file), verbose = F, 
                           clinicalData = annot)
    maf_lung = merge_mafs(list(maf_lung, current_maf),verbose = F,cnTable = gl_df_lung)
  }
}

# C. MAF summary 
# head(maf_lung@data, 2) # MAF similar to a VCF
getSampleSummary(maf_lung)

variant_cols = c("#33A02C","#1f78b4","#000000","#535c68","#efa9ae","#f4e87c","indianred1")
names(variant_cols) = c("Missense_Mutation","Frame_Shift_Del","Multi_Hit","pathway","GL missense",      
                        "GL splice acceptor","In_Frame_Del")

plotmafSummary(maf = maf_lung, rmOutlier = F, addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE, color = variant_cols,showBarcodes = T)

maf_onco = subsetMaf(maf_lung, query = "!grepl('Blood', Tumor_Sample_Barcode)")

if(plot_figs) png(file.path(resdir,"MAF_lung_noannot.png"),  width = 1000,height = 750,res = 300)
oncoplot(maf = maf_onco, showTumorSampleBarcodes = F, sampleOrder =  annot$Sample[annot$Patient=="Lung"], 
         cohortSize = 4,clinicalFeatures = "Amount", annotationFontSize = 0.8, fontSize = 0.6,
         legendFontSize = 0, legend_height = 2, gene_mar = 4, barcode_mar = 6, drawRowBar = F,
         writeMatrix = F, drawBox = F, anno_height = .8, titleText = "", colors = variant_cols
)
if(plot_figs) dev.off()

if(plot_figs) png(file.path(resdir,"MAF_lung_withannot.png"),  width = 10000,height = 7500,res = 1500)
oncoplot(maf = maf_onco, showTumorSampleBarcodes = F, sampleOrder =  annot$Sample[annot$Patient=="Lung"], 
         cohortSize = 4,clinicalFeatures = "Amount", annotationFontSize = 0.5, fontSize = 0.6,
         legendFontSize = 0.5, legend_height = 2, gene_mar = 4, barcode_mar = 6, drawRowBar = F,
         writeMatrix = F, drawBox = F, anno_height = .8, titleText = "", colors = variant_cols
)
if(plot_figs) dev.off()

# none of these genes are seen as top lung cancer genes in literature

# 10. Plasma vs Urine ----

IDs_list = c() # generate list of mutation IDs present in each sample 
for(i in 1:length(unique(vcf_somatic$Sample))){
  samp = unique(vcf_somatic$Sample)[i]
  current_IDs = vcf_somatic %>% 
    filter(Sample == samp) %>% 
    pull(ID)
  
  IDs_list = c(IDs_list, list(current_IDs))
  names(IDs_list)[i] = samp 
}

Patient_1_urine = unlist(IDs_list[annot$Sample[annot$Patient == "Patient_1" & annot$Tissue == "Urine"]])
Patient_1_plasma = unlist(IDs_list[annot$Sample[annot$Patient == "Patient_1" & annot$Tissue == "Plasma"]])

if(plot_figs) png(file.path(resdir,"Patient_1_venn.png"), width = 1000,height = 1000,res = 300)
venn(x = list("Urine" = Patient_1_urine, "Plasma" = Patient_1_plasma),
     zcolor = "style", box = F, ggplot = T)+
  theme(plot.background = element_blank(), panel.background = element_blank())
if(plot_figs) dev.off()

# no. of unique somatic mutations in Patient 1
length(unique(vcf_somatic$ID[grepl("1",vcf_somatic$Patient)]))

# consequence of Patient exonic mutations only seen in urine  
urine_only = vcf_somatic %>% 
  filter(!ID %in% Patient_1_plasma & ID %in% Patient_1_urine, 
         Patient == "Patient_1", Tissue == "U") %>% 
  select(ID, Gene.refGene, ExonicFunc.refGene) %>% 
  distinct() 

urine_only %>% 
  count(ExonicFunc.refGene)

# consequence of Patient exonic mutations only seen in plasma   
plasma_only = vcf_somatic %>% 
  filter(ID %in% Patient_1_plasma & !ID %in% Patient_1_urine, 
         Patient == "Patient_1", Tissue == "PLASMA") %>% 
  select(ID, Gene.refGene, ExonicFunc.refGene) %>% 
  distinct() 

plasma_only %>% 
  count(ExonicFunc.refGene)

plasma_only %>% 
  filter(Gene.refGene %in% cancer_genes)

urine_only %>% 
  filter(Gene.refGene %in% cancer_genes)

# oncoplot consistency 
p1_oncofreq = maf_combined@data %>% 
  as.data.frame() %>% 
  select(where(~!all(is.na(.x)))) %>% 
  select(Tumor_Sample_Barcode, Gene = Hugo_Symbol) %>% 
  left_join(annot) %>% 
  filter(Patient=="Patient_1") %>% 
  select(Sample, Patient, Timepoint, Tissue, Gene) %>% 
  count(Tissue, Gene) %>% 
  filter(Tissue != "Blood", Gene %in% oncoplot_genes) %>% 
  arrange(Gene,Tissue) %>%
  mutate(n_samps = ifelse(Tissue=="Plasma",5,4),
         prop = n/n_samps,
         prop = ifelse(prop > 1, 1, prop),
         mostmuts = NA) %>% 
  group_by(Gene)

for(mygene in unique(p1_oncofreq$Gene)){
  myprops = p1_oncofreq$prop[p1_oncofreq$Gene == mygene]
  if(length(myprops)==1){
    p1_oncofreq$mostmuts[p1_oncofreq$Gene == mygene] = NA
  }else{
    if(myprops[1] > myprops[2]){
      p1_oncofreq$mostmuts[p1_oncofreq$Gene == mygene] = "Plasma"
    }
    if(myprops[1] < myprops[2]){
      p1_oncofreq$mostmuts[p1_oncofreq$Gene == mygene] = "Urine"
    }
    if(myprops[1] == myprops[2]){
      p1_oncofreq$mostmuts[p1_oncofreq$Gene == mygene] = "Equal"
    }
  }
}
table(p1_oncofreq$mostmuts) / 2 # whether a gene in the oncoplot is seen more frequently in plasma, urine or equal


# TMB 
if(plot_figs) png(file.path(resdir,"TMB_patient1.png"), width = 1300,height = 1000,res = 300)
tmb_df %>% 
  filter(Patient == "Patient_1") %>% 
  group_by(Timepoint,Tissue) %>% 
  mutate(tpointmean = mean(TMB)) %>% 
  ggplot() + 
    geom_line( aes(x = n_weeks, y = tpointmean, group = Tissue)) + 
    geom_point(aes(x = n_weeks, y = tpointmean, colour = Tissue), size = 2.5) + 
    scale_x_continuous(breaks = c(0,6,12,18,24,30), labels = c("Diagnosis",6,12,18,24,30)) + 
    scale_y_continuous(limits = c(0,27), expand = c(0,0)) +
    theme_bw() + theme(legend.title = element_markdown()) + 
    scale_color_manual(drop = F, values = tissue_cols, name = "**Patient 1**") + 
    labs(y = "Tumour mutation burden (mut/Mb)", x = "Weeks after diagnosis")
dev.off()

tmb_df %>% 
  filter(Patient == "Patient_1") %>% 
  group_by(Tissue, Timepoint) %>% 
  mutate(tpointmean = mean(TMB)) %>% 
  ungroup() %>% 
  group_by(Timepoint) %>%  
  mutate(minval = min(tpointmean),
         maxval = max(tpointmean),
         perdiff = (maxval - minval) * 100 / ((minval + maxval)/ 2)
         ) %>% 
  arrange(Timepoint) %>% 
  select(Timepoint,maxval, minval, perdiff) %>% 
  distinct()


# 11. Mutational Signatures ----
library(Palimpsest) # devtools::install_github("FunGeST/Palimpsest")
library(BSgenome.Hsapiens.UCSC.hg19) # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

plot_sigs_edit = function (input_data = NULL, Title = NA, label = "Full") {
  requireNamespace("scales", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  if (label %!in% c("Full", "Top", "Bottom", "None")) 
    stop("label must be 'Full', 'Top', 'Bottom' or 'None'")
  Individual <- ifelse(1 %in% dim(input_data), TRUE, FALSE)
  if (Individual == FALSE) 
    num_of_sigs = nrow(input_data)
  if (Individual == TRUE) 
    num_of_sigs = 1
  input_data <- as.matrix(input_data)
  for (i in 1:num_of_sigs) {
    if (Individual == FALSE) {
      if (ncol(input_data) %!in% c(78, 83, 96)) 
        stop("input_data format incorrect")
      if (ncol(input_data) == 96) 
        Type <- "SBS"
      if (ncol(input_data) == 78) 
        Type <- "DBS"
      if (ncol(input_data) == 83) 
        Type <- "ID"
      Context <- colnames(input_data)
      plot_input <- as.data.frame(input_data[i, ])
      if (length(Title) == 1) 
        plot_title <- paste(Title, i, sep = ".")
      if (length(Title) > 1) 
        plot_title <- Title[i]
      if (is.na(Title[1])) 
        plot_title <- ""
    }
    if (Individual == TRUE) {
      if (nrow(input_data) == 1) 
        input_data <- t(input_data)
      if (nrow(input_data) %!in% c(38, 78, 83, 96)) 
        stop("input_data format incorrect")
      if (nrow(input_data) == 96) 
        Type <- "SBS"
      if (nrow(input_data) == 78) 
        Type <- "DBS"
      if (nrow(input_data) == 83) 
        Type <- "ID"
      if (nrow(input_data) == 38) 
        Type <- "SV"
      Context <- rownames(input_data)
      plot_input <- data.frame(input_data)
      if (!is.na(Title)) 
        plot_title <- Title
      if (is.na(Title)) 
        plot_title <- ""
    }
    colnames(plot_input)[1] <- "freq"
    max.y = max(plot_input$freq)
    
    if (Type == "SBS") {
      plot_input = plot_input %>% mutate(Substype = paste0(substr(Context, 
                                                                  1, 1), ">", substr(Context, 2, 2)), Context = paste0(substr(Context, 
                                                                                                                              4, 4), substr(Context, 1, 1), substr(Context, 
                                                                                                                                                                   6, 6)), Substype_blank = NA, Colours = c(rep("skyblue3", 
                                                                                                                                                                                                                16), rep("black", 16), rep("red", 16), rep("grey", 
                                                                                                                                                                                                                                                           16), rep("green", 16), rep("pink", 16)), Context_numerical = c(1:nrow(plot_input)), 
                                         Start_pos = c(rep(1, 16), rep(17, 16), rep(33, 
                                                                                    16), rep(49, 16), rep(65, 16), rep(81, 16)), 
                                         Substype_length = 16, Substype_numerical = c(as.numeric(as.factor(Substype))))
      prev <- "rien de rien"
      for (j in 1:nrow(plot_input)) {
        if (plot_input$Substype[j] != prev) {
          plot_input$Substype_blank[j] <- plot_input$Substype[j]
          prev <- plot_input$Substype[j]
          next
        }
        else {
          plot_input$Substype_blank[j] <- paste(" ")
        }
      }
      xmax_correction = 20
    }
    
    max.y = max(plot_input$freq)
    plot_output = ggplot(data = plot_input) + geom_bar(aes(x = as.factor(plot_input$Context_numerical), 
                                                           y = freq), stat = "identity", fill = plot_input$Colours, 
                                                       width = 0.5) + ggtitle(plot_title) + ylab("") + 
      xlab("") + coord_cartesian(xlim = c(1, nrow(plot_input)), 
                                 ylim = c(0, max.y * 1.3), clip = "off") + theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), 
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_blank(), panel.background = element_blank(), 
            plot.background = element_blank(), axis.text.y = element_text(size = 15, 
                                                                          colour = "black"), axis.ticks.x = element_blank(), 
            plot.margin = unit(c(1.2, 1.6, 1.2, 1.1), "cm"), 
            plot.title = element_text(size = 55, face = "bold", 
                                      hjust = 0.02, vjust = -5))
    
    if (Type == "SBS") {
      if (label == "Full") {
        max.y_lab <- max.y
        if ((round(max.y_lab, 2) * 100)%%2 == 1) 
          max.y_lab <- round(max.y_lab, 2) + 0.01
        plot_output = plot_output + xlab(" ") + ylab(" ") + 
          scale_x_discrete("", labels = plot_input$Context) + 
          annotate("text", x = -3.5, y = max.y * 0.5, 
                   label = "Proportion of Mutations", angle = 90, 
                   size = 7, fontface = "plain", colour = "black") + 
          annotate("text", x = 48, y = 0 - max.y * 0.22, 
                   label = "Trinucleotide Context", size = 8, 
                   fontface = "plain", colour = "black") + 
          theme(axis.text.x = element_text(angle = 90, 
                                           vjust = 0.5, hjust = 0, size = 20, colour = "black", 
                                           family = "mono"), axis.text.y = element_text(size = 20), 
                axis.ticks.y.right = element_blank()) + 
          geom_rect(mapping = aes(xmin = (plot_input$Start_pos - 
                                            0.5), ymin = max.y * 1.25, xmax = (plot_input$Start_pos + 
                                                                                 xmax_correction), ymax = max.y * 1.375), 
                    fill = plot_input$Colours, color = "white", 
                    size = 2.5) + geom_text(position = "identity", 
                                            mapping = aes(x = plot_input$Start_pos + (plot_input$Substype_length/2.2), 
                                                          y = max.y * 1.42, label = plot_input$Substype_blank), 
                                            size = 10) + scale_y_continuous(expand = c(0, 
                                                                                       0), limits = c(-1, max.y * 1.42), breaks = seq(from = 0, 
                                                                                                                                      to = max.y_lab, by = max.y_lab/2), labels = scales::number_format(accuracy = 0.01), 
                                                                            sec.axis = sec_axis(~. + 0, labels = NULL)) + 
          geom_rect(mapping = aes(xmin = (96.5), ymin = max.y * 
                                    1.25, xmax = 100, ymax = max.y * 1.375), 
                    fill = "white", color = "white", size = 2.5)
      }
    }
    
    plot(plot_output)
  }
}

vcf_palimp = annotate_VCF(vcf = vcf_somatic)

SBS_input = palimpsest_input(vcf = vcf_palimp,Type ="SBS")

inputsamples = annot$Sample[grepl("1",annot$Patient)]

plotin = SBS_input$mut_nums[,colnames(SBS_input$mut_nums) %in% inputsamples] %>% 
  rowSums() %>% 
  as.matrix()

plotin = plotin / sum(plotin)

if(plot_figs) png(file.path(resdir,"96_mutcat_Patient_1.png"), width = 2400,height = 700,res = 100)
plot_sigs_edit(plotin, Title = "Patient 1")
if(plot_figs) dev.off()

patient1_cos_sim = compare_results(reference_sigs = SBS_cosmic,extraction_1 = t(plotin), extraction_1_name = "Patient_1",
                lower_threshold = 0) %>% 
  arrange(desc(Ref_Patient_1_cosine_score)) %>% 
  mutate(sigrank = c(1:nrow(.)))
# most similar to signature 55: https://cancer.sanger.ac.uk/signatures/sbs/sbs55/

MSI_sigs = paste0("SBS",c(6, 15, 20, 21, 26))

report_cosine_tab = head(patient1_cos_sim, 3) %>% 
  rbind(filter(patient1_cos_sim, Ref_Signature %in% MSI_sigs)) %>% 
  select(`COSMIC Signature` = Ref_Signature, `Cosine similiartity` = Ref_Patient_1_cosine_score,
         Rank = sigrank)

write_delim(report_cosine_tab, 
            file = file.path(resdir,("Patient_1_COSMIC_cosines.txt")),
            delim = "\t",col_names = T)

