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
#           User Input
###########################################################
## Phenotype file
phenotype <- read.csv("/path/to/cognitive_data.csv")
## (sparse) GRM file
sgrm <- get(load("/path/to/Q0_unre_Caucasian_chrall_pruned.sparseGRM.sGRM.RData"))
sgrm@Dimnames[[1]]<-substr(sgrm@Dimnames[[1]], 3, 9)
sgrm@Dimnames[[2]]<-substr(sgrm@Dimnames[[2]], 3, 9)
sgrm@Dimnames[[1]]<-paste(sgrm@Dimnames[[1]],sgrm@Dimnames[[1]],sep = "_")
sgrm@Dimnames[[2]]<-paste(sgrm@Dimnames[[2]],sgrm@Dimnames[[2]],sep = "_")

## file directory for the output file 
output_path <- "/path/to/STEP_1/"
## output file name
output_name <- "obj_nullmodel.Rdata"

###########################################################
#           Main Function 
###########################################################
## fit null model
obj_nullmodel <- fit_nullmodel(pheno~age+as.factor(sex)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,
                               data=phenotype,kins=sgrm,use_sparse=TRUE,kins_cutoff=0.088,id="IID",use_SPA=FALSE,
                               family=gaussian(link="identity"),verbose=TRUE)

save(obj_nullmodel,file=paste0(output_path,output_name))
