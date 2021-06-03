# Genetic gain19 germination tests
#combine raw data into a single file and standardize phenotypes

library(readxl)
library(tidyverse)
library(readxl)

#this is where data has been stored; could change this? 
setwd(dir = rprojroot::find_rstudio_root_file())

file_list <- list.files(path='PhenotypeData/RawData/2019/')
TP1=which(substr(file_list,8,10)=="TP1")
TP2=which(substr(file_list,8,10)=="TP2")
TP3=which(substr(file_list,8,10)=="TP3")
file_list = paste0('PhenotypeData/RawData/2019/',file_list)

#quickly and efficiently merge data by timepoint. Data could be merged into a single file as well. 
TP1all <- data.frame()
for (i in TP1){
  temp_data <- read_excel(file_list[i], range = cell_cols("A:N")) 
  temp_data$TP <- "TP1" 
  TP1all <- rbind(TP1all, temp_data) 
}

TP2all <- data.frame()
for (i in TP2){
  temp_data <- read_excel(file_list[i], range = cell_cols("A:O")) 
  temp_data$TP <- "TP2" 
  TP2all <- rbind(TP2all, temp_data) 
}

TP3all_30k <- data.frame()
for (i in TP3){
  temp_data <- read_excel(file_list[i],sheet=1, range = cell_cols("A:O")) 
  temp_data$TP <- "TP3" 
  TP3all_30k <- rbind(TP3all_30k, temp_data) 
}

TP3all_100k <- data.frame()
for (i in TP3){
  temp_data <- read_excel(file_list[i],sheet=2, range = cell_cols("A:N")) 
  temp_data$TP <- "TP3" 
  TP3all_100k <- rbind(TP3all_100k, temp_data) 
}

TP1all[,11:15] = NA #remove WS data for TP1
colnames(TP1all)[11:15]=colnames(TP2all)[11:15]
TP1all$TP="TP1"

#stack the two reps into a single column for each day
TP1all=rbind(TP1all,TP1all)
TP1all[665:1328,6:10]=TP1all[665:1328,11:15]
TP1all=TP1all[,-c(11:15)]; TP1all$rep=rep(c("A","B"), each=664)
TP1all$TP="TP1"

TP2all=rbind(TP2all,TP2all)
TP2all[665:1328,6:10]=TP2all[665:1328,11:15]
TP2all=TP2all[,-c(11:15)]; TP2all$rep=rep(c("A","B"), each=664)

TP3all_30k=rbind(TP3all_30k,TP3all_30k)
TP3all_30k[665:1328,6:10]=TP3all_30k[665:1328,11:15]
TP3all_30k=TP3all_30k[,-c(11:15)]; TP3all_30k$rep=rep(c("A","B"), each=664)


#convert data to traits

Germ_traits=function(df){
    df$GE=apply(df[,6:9],1,function(x) {sum(as.numeric(x[1:3]))/sum(as.numeric(x[1:4]))})  #does not count dead kernels in total count
    df$GC=apply(df[,6:10],1,function(x) {sum(as.numeric(x[1:4]))/sum(as.numeric(x[1:5]))}) #does count dead kernels
    df$GI=apply(df[,6:8],1, function(x) {(sum(as.numeric(x[1:3]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3])))*10})
    df$GIscale=apply(df[,c(13,15)],1, function(x) {as.numeric(x[1])*as.numeric(x[2])}) #GI*GE
    df$StartDate=as.factor(df$StartDate)
    df=df[-c(which(df$Entry=="JIC")),]
    df$rep=as.factor(df$rep)
    df$TP=as.factor(df$TP)
    df
}

Germ_traits_100k=function(df){
  df$GE=apply(df[,6:9],1,function(x) {sum(as.numeric(x[1:3]))/sum(as.numeric(x[1:4]))})  #does not count dead kernels in total count
  df$GC=apply(df[,6:10],1,function(x) {sum(as.numeric(x[1:4]))/sum(as.numeric(x[1:5]))}) #does count dead kernels
  df$GI=apply(df[,6:8],1, function(x) {(sum(as.numeric(x[1:3]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3])))*10})
  df$GIscale=apply(df[,c(16,18)],1, function(x) {as.numeric(x[1])*as.numeric(x[2])}) #GI*GE
  df$WS=apply(df[,c(6:8,11:13)], 1, function(x) {sum(as.numeric(x[1:3])) - sum(as.numeric(x[4:6]))})
  df$WSI=apply(df[,11:13],1, function(x) {(sum(as.numeric(x[1:3]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3])))*10})
  df$StartDate=as.factor(df$StartDate)
  df=df[-c(which(df$Entry=="JIC")),]
  df$TP=as.factor(df$TP)
  df
}

TP1all=Germ_traits(TP1all)
TP2all=Germ_traits(TP2all)
TP3all_30k=Germ_traits(TP3all_30k) 
TP3all_100k=Germ_traits_100k(TP3all_100k) 

TP1all[which(TP1all$GI=='NaN'),15:16]=NA

#For some reason Pinnacle was spelled wrong; needs to match the spelling in genotype file
TP2all$Entry[which(TP2all$Entry=="Pincle")]="Pinnacle"
TP3all_30k$Entry[which(TP3all_30k$Entry=="Pincle")]="Pinnacle"
TP3all_100k$Entry[which(TP3all_100k$Entry=="Pincle")]="Pinnacle"

# save(TP1all, file="TP1 all data for analysis.RData")
# save(TP2all, file="TP2 all data for analysis.RData")
# save(TP3all_30k, file="TP3 30k all data for analysis.RData")
# save(TP3all_100k, file="TP3 100k all data for analysis.RData")