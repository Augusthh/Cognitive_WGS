##################################################################################
# Gene-centric analysis for coding rare variants in long masks using STAARpipeline
# Xihao Li, Zilin Li
##################################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)


###########################################################
#           Project paths
###########################################################
project_root <- Sys.getenv("COGNITIVE_WGS_ROOT")
data_dir <- file.path(project_root, "data")
staar_dir <- file.path(project_root, "STAARpipeline")
result_dir <- file.path(project_root, "result")

###########################################################
#           User Input
###########################################################
pheno_name <- commandArgs(TRUE)[2]
## aGDS directory
agds_dir <- get(load(file.path(staar_dir, "STEP_0", "agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(file.path(result_dir, pheno_name, "STEP_1", "obj_nullmodel.Rdata")))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load(file.path(staar_dir, "STEP_0", "Annotation_name_catalog.Rdata")))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- file.path(result_dir, pheno_name, "STEP_3.1")
## output file name
output_file_name <- paste0(pheno_name, "_coding")
## input array id from batch file (Harvard FAS RC cluster)
arrayid_longmask <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

## analyze large coding masks
arrayid <- c(57,112,112,113,113,113,113,113,113,113)
sub_seq_id <- c(840,543,544,575,576,577,578,579,580,582)

region_spec <- data.frame(arrayid,sub_seq_id) 
sub_seq_id <- ((arrayid_longmask-1)*5+1):min(arrayid_longmask*5,length(arrayid))

### aGDS file
genes <- genes_info

results_coding <- c()
for(kk in sub_seq_id)
{
  print(kk)
  arrayid <- region_spec$arrayid[kk]
  sub_id <- region_spec$sub_seq_id[kk]
  
  chr <- which.max(arrayid <= cumsum(group.num.allchr))
  agds.path <- agds_dir[chr]
  genofile <- seqOpen(agds.path)
  
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  gene_name <- genes_info_chr[sub_id,1]
  
  results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_coding <- append(results_coding,results)
  
  seqClose(genofile)
}

save(results_coding, file = file.path(output_path, paste0(output_file_name, "_", arrayid_longmask + 379, ".Rdata")))

