###########################################################
# fit STAAR null model
# Xihao Li, Zilin Li
###########################################################
rm(list=ls())
gc()

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
## Phenotype file
pheno_name <- commandArgs(TRUE)[1]
phenotype <- read.csv(file.path(data_dir, "cognitive_10phenotype.csv"))
## (sparse) GRM file
sgrm <- get(load(file.path(staar_dir, "Q0_unre_Caucasian_chrall_pruned.sparseGRM.sGRM.RData")))
sgrm@Dimnames[[1]]<-substr(sgrm@Dimnames[[1]], 3, 9)
sgrm@Dimnames[[2]]<-substr(sgrm@Dimnames[[2]], 3, 9)
sgrm@Dimnames[[1]]<-paste(sgrm@Dimnames[[1]],sgrm@Dimnames[[1]],sep = "_")
sgrm@Dimnames[[2]]<-paste(sgrm@Dimnames[[2]],sgrm@Dimnames[[2]],sep = "_")

## file directory for the output file 
output_path <- file.path(result_dir, pheno_name, "STEP_1")
if (!dir.exists(output_path)) {dir.create(output_path, recursive = TRUE)}
## output file name
output_name <- "obj_nullmodel.Rdata"

###########################################################
#           Main Function 
###########################################################
## fit null model
trait_specific_cov <- commandArgs(TRUE)[2]
formula=paste0(pheno_name,'~age+as.factor(sex)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',trait_specific_cov)
obj_nullmodel <- fit_nullmodel(formula,
                               data=phenotype,kins=sgrm,use_sparse=TRUE,kins_cutoff=0.088,id="IID",use_SPA=FALSE,
                               family=gaussian(link="identity"),verbose=TRUE)

save(obj_nullmodel, file = file.path(output_path, output_name))
