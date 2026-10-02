###########################################################
# Individual analysis using STAARpipeline
# Xihao Li, Zilin Li
###########################################################
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
## Number of jobs for each chromosome
jobs_num <- get(load(file.path(staar_dir, "STEP_0", "jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(file.path(staar_dir, "STEP_0", "agds_dir.Rdata")))
## Null model
pheno_name <- commandArgs(TRUE)[2]
obj_nullmodel <- get(load(file.path(result_dir, pheno_name, "STEP_1", "obj_nullmodel.Rdata")))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
output_path <- file.path(result_dir, pheno_name, "STEP_2")
if (!dir.exists(output_path)) {dir.create(output_path, recursive = TRUE)}
## output file name
output_file_name <- paste0(pheno_name,"_Individual_Analysis")
## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis, file = file.path(output_path, paste0(output_file_name, "_", arrayid, ".Rdata")))

seqClose(genofile)
