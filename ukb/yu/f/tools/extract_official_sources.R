#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(readxl); library(data.table) })
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop("Usage: Rscript extract_official_sources.R <supplement.xlsx> <project_dir>")
x<-normalizePath(args[[1]],mustWork=TRUE); root<-normalizePath(args[[2]],mustWork=TRUE)
d<-as.data.table(read_excel(x,sheet="S12",skip=1))
cad<-d[Disease=="CAD" & !is.na(Proteins),.(protein=toupper(Proteins),rank=as.integer(Rank),published_importance=as.numeric(`Protein importance for one disease`))]
s23<-as.data.table(read_excel(x,sheet="S23",skip=1))[,.(protein=toupper(Protein),uniprot=UniProt,panel=Panel,protein_definition=Protein_definition)]
cad<-merge(cad,s23,by="protein",all.x=TRUE,sort=FALSE); setorder(cad,rank)
cad[,local_feature:=gsub("^_|_$","",gsub("[^a-z0-9]+","_",tolower(protein)))]
stopifnot(nrow(cad)==257L,uniqueN(cad$protein)==257L,!anyNA(cad$uniprot))
fwrite(cad,file.path(root,"config","yu_cad_257_official.csv"))
s10<-as.data.table(read_excel(x,sheet="S10",skip=1)); i<-which(s10$Disease=="CAD")
ref<-s10[(i+1):(i+3),.(model=Model,auc_published=AUC,accuracy_published=Accuracy,sensitivity_published=Sensitivity,specificity_published=Specificity,f1_published=F1.score,brier_published=BRIER.score)]
fwrite(ref,file.path(root,"config","yu_cad_published_metrics.csv"))
cat("Extracted",nrow(cad),"official CAD proteins and",nrow(ref),"reference models.\n")
