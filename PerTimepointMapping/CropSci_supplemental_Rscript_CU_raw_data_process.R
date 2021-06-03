#### Supplemental R script for raw CU/JIC seed germination data processing for "Pleiotropic effects of HvMKK3 alleles on preharvest sprouting, seed dormancy, and malting quality"
setwd(rprojroot::find_rstudio_root_file())
getwd()


#load geno data
load("GenotypeData/myGD_LDpruned_w_KASP.RData")
load("GenotypeData/myGM_LDpruned_w_KASP.RData")
#load 2020 raw pheno data
load("PhenotypeData/RawData/2020/GGS_TP1_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP2_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP3_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP4_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP5_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP6_for_analysis.RData")
load("PhenotypeData/RawData/2020/GGS_TP7_for_analysis.RData")


#preliminary edits
GGS_TP3$Day2Germ=as.numeric(GGS_TP3$Day2Germ)

germ_phenos = function(df){
  df$GE3=round(apply(df[,8:13],1,function(x){sum(x[1:3])/sum(x[1:6])}),2)
  df$GE5=round(apply(df[,8:13],1,function(x){sum(x[1:5])/sum(x[1:6])}),2)
  
  df$GI3=round(apply(df[,8:13],1,function(x){(sum(as.numeric(x[1:3]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3])))*10}),2)
  df$GI5=round(apply(df[,8:13],1,function(x){(sum(as.numeric(x[1:5]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3]) + 4*as.numeric(x[4]) + 5*as.numeric(x[5])))*10}),2)
  df$GI3[which(apply(df[,8:10],1,sum)==0)] = 0
  
  df$GI3scale=df$GI3*df$GE3
  df$GI5scale=df$GI5*df$GE5
  df$rep=as.factor(df$rep)
  df$Location=as.factor(df$Location)
  df$CoolerLetter=as.factor(df$CoolerLetter)
  
  #standardize Entry names with geno file
  df$Entry=gsub("\\.","_",df$Entry)
  df$Entry=gsub("\\-","_",df$Entry)
  df$Entry=gsub(" ","_",df$Entry)
  
  
  df
} #create seed germination traits and process raw data

GGS_TP1=germ_phenos(GGS_TP1)
GGS_TP2=germ_phenos(GGS_TP2)
GGS_TP3=germ_phenos(GGS_TP3)
GGS_TP4=germ_phenos(GGS_TP4)
GGS_TP5=germ_phenos(GGS_TP5)
GGS_TP6=germ_phenos(GGS_TP6)

#address outliers 
GGS_TP1[152,16:21] = NA #incomplete data entry
GGS_TP5[c(2204,904),16:21] = NA #very low germ, likely due to disease
GGS_TP6[which(GGS_TP6$KernelsLeft>5),16:21] =NA #very low germ, likely due to disease


##### load 2019 raw phenotype data #####
#IMPORTANT: 2019 TP1 matches with 2020 TP1, 2019 TP2 matches with 2020 TP4, and 2019 TP3 matches with 2020 TP6
load("PhenotypeData/RawData/2019/TP1 all data for analysis.RData")
load("PhenotypeData/RawData/2019/TP2 all data for analysis.RData")
load("PhenotypeData/RawData/2019/TP3 30k all data for analysis.RData");TP3all=TP3all_30k; rm(TP3all_30k)


##### heritability estimates #####
library(lme4)
Cullis_H2=function(model){
  library(arm)
  ses<- se.ranef(model)$'Entry' #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(model, comp="Variance")$'Entry'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated with the model saved in 'm'
  H2<- round(mean(Reliability),3) #This is equivalent to broad-sense heritability on the line-mean (or family-mean, if your individuals are non-inbred families) basis
  H2
}
  #TP1
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),])) #0.969
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),])) #0.983
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),])) #0.988
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),])) #0.989
  #TP2
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),])) #0.96
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),])) #0.976
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),])) #0.983
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),])) #0.984
  #TP3
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),])) #0.947
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),])) #0.951
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),])) #0.97
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),])) #0.97
  #TP4
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),])) #0.927
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),])) #0.949
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),])) #0.965
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),])) #0.964
  #TP5
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),])) #0.674
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),])) #0.779
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),])) #0.917
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),])) #0.918
  #TP6
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),])) # .687 
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),])) #.69
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),])) #.895
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),])) #.896
  

##### fixed effects models #####
tp1_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif; row 152 has incomplete data and is an outlier
tp1_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif
tp1_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif
tp1_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif

tp2_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif
tp2_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif

tp3_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #rep not signif
tp3_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #rep not signif

tp4_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #rep not signif
tp4_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #rep not signif

tp5_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif
tp5_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif

tp6_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif
tp6_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif


##### random effects models #####
library(lme4)
tp1_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif; row 152 has incomplete data and is an outlier
tp1_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif
tp1_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif
tp1_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif

tp2_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif
tp2_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif

tp3_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #rep not signif
tp3_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #rep not signif

tp4_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #rep not signif
tp4_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #rep not signif

tp5_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif
tp5_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif

tp6_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif
tp6_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif


##### BLUEs for each trait, organized by timepoint #####
get_entry_fixed_Effects = function(model, trait){
  x = data.frame(coef(model))
  x$names = rownames(x)
  x$names2 = substr(x$names,1,5)
  x = subset(x, x$names2 %in% c('(Inte', 'Entry'))
  x$taxa = substring(x$names, 6)
  x$taxa[1] = model$xlevels$Entry[1]
  intercept = x$coef.model.[1]
  x$coef.model.[1] = 0
  x$coef.model. = x$coef.model.+intercept
  x = data.frame(taxa = x$taxa, x$coef.model.)
  colnames(x)[2] = trait
  return(x)
} #need this to run 'BLUE_summary'
BLUE_summary=function(mods,PM,TP,date){
  l=list(get_entry_fixed_Effects(mods[[1]], 'GE3'),get_entry_fixed_Effects(mods[[2]], 'GE5'),
         get_entry_fixed_Effects(mods[[3]], 'GI3scale'),get_entry_fixed_Effects(mods[[4]], 'GI5scale'))
  
  TP_summary <- Reduce(
    function(x, y, ...) merge(x, y, by="taxa" ,all = TRUE),
    l
  )
  TP_summary=TP_summary[-c(which(!TP_summary$taxa %in% myGD20_prune$taxa)),]
  TP_summary$PM_date=PM
  TP_summary$TP=TP
  TP_summary$date=date
  TP_summary
} #combine trait BLUEs across timepoints

TP1_GGS20_BLUE_all=BLUE_summary(list(tp1_ge3.lm,tp1_ge5.lm,tp1_gi3.lm,tp1_gi5.lm), PM=5, TP="TP1", date="9/27/2020")
TP2_GGS20_BLUE_all=BLUE_summary(list(tp2_ge3.lm,tp2_ge5.lm,tp2_gi3.lm,tp2_gi5.lm), PM=19, TP="TP2", date="10/19/2020")
TP3_GGS20_BLUE_all=BLUE_summary(list(tp3_ge3.lm,tp3_ge5.lm,tp3_gi3.lm,tp3_gi5.lm), PM=33, TP="TP3", date="10/26/2020")
TP4_GGS20_BLUE_all=BLUE_summary(list(tp4_ge3.lm,tp4_ge5.lm,tp4_gi3.lm,tp4_gi5.lm), PM=47, TP="TP4", date="11/9/2020")
TP5_GGS20_BLUE_all=BLUE_summary(list(tp5_ge3.lm,tp5_ge5.lm,tp5_gi3.lm,tp5_gi5.lm), PM=68, TP="TP5", date="11/30/2020")
TP6_GGS20_BLUE_all=BLUE_summary(list(tp6_ge3.lm,tp6_ge5.lm,tp6_gi3.lm,tp6_gi5.lm), PM=110, TP="TP6", date="1/11/2021")


##### BLUPs for each trait, organized by timepoint #####
gapit_blup=function(mod,trait){
  df=ranef(mod)$'Entry' + fixef(mod)[1]
  df=cbind(rownames(df), df); 
  colnames(df)=c("taxa", trait); 
  df[,1]=as.character(df[,1])
  df
} #need this to run 'BLUP_summary'
BLUP_summary=function(mods,PM,TP,date){
  l=list(gapit_blup(mods[[1]], 'GE3'),gapit_blup(mods[[2]], 'GE5'),
         gapit_blup(mods[[3]], 'GI3scale'),gapit_blup(mods[[4]], 'GI5scale'))
  
  TP_summary <- Reduce(
    function(x, y, ...) merge(x, y, by="taxa" ,all = TRUE),
    l
  )
  TP_summary=TP_summary[-c(which(!TP_summary$taxa %in% myGD20_prune$taxa)),]
  TP_summary$PM_date=PM
  TP_summary$TP=TP
  TP_summary$date=date
  TP_summary
} #combine trait BLUPs across timepoints
TP1_GGS20_BLUP_all=BLUP_summary(list(tp1_ge3.lmer,tp1_ge5.lmer,tp1_gi3.lmer,tp1_gi5.lmer), PM=5, TP="TP1", date="9/27/2020")
TP2_GGS20_BLUP_all=BLUP_summary(list(tp2_ge3.lmer,tp2_ge5.lmer,tp2_gi3.lmer,tp2_gi5.lmer), PM=19, TP="TP2", date="10/19/2020")
TP3_GGS20_BLUP_all=BLUP_summary(list(tp3_ge3.lmer,tp3_ge5.lmer,tp3_gi3.lmer,tp3_gi5.lmer), PM=33, TP="TP3", date="10/26/2020")
TP4_GGS20_BLUP_all=BLUP_summary(list(tp4_ge3.lmer,tp4_ge5.lmer,tp4_gi3.lmer,tp4_gi5.lmer), PM=47, TP="TP4", date="11/9/2020")
TP5_GGS20_BLUP_all=BLUP_summary(list(tp5_ge3.lmer,tp5_ge5.lmer,tp5_gi3.lmer,tp5_gi5.lmer), PM=68, TP="TP5", date="11/30/2020")
TP6_GGS20_BLUP_all=BLUP_summary(list(tp6_ge3.lmer,tp6_ge5.lmer,tp6_gi3.lmer,tp6_gi5.lmer), PM=110, TP="TP6", date="1/11/2021")

##### Combined 2019/2020 analysis #####
####  Only for GE3 and GI3
TP1all$Entry=gsub("\\-","_",TP1all$Entry)
TP1all$Entry=gsub(" ","_",TP1all$Entry)
TP1all$Entry=gsub("\\.","_",TP1all$Entry)

TP2all$Entry=gsub("\\-","_",TP2all$Entry)
TP2all$Entry=gsub(" ","_",TP2all$Entry)
TP2all$Entry=gsub("\\.","_",TP2all$Entry)

TP3all$Entry=gsub("\\-","_",TP3all$Entry)
TP3all$Entry=gsub(" ","_",TP3all$Entry)
TP3all$Entry=gsub("\\.","_",TP3all$Entry)

TP1all$Env=ifelse(substr(TP1all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP2all$Env=ifelse(substr(TP2all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP3all$Env=ifelse(substr(TP3all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP1all$rep=ifelse(TP1all$rep=="A", "1","2")
TP2all$rep=ifelse(TP2all$rep=="A", "1","2")
TP3all$rep=ifelse(TP3all$rep=="A", "1","2")

colnames(TP1all)[c(4,13,16)]=c("Plot", "GE3","GI3scale")
colnames(TP2all)[c(4,13,16)]=c("Plot","GE3","GI3scale")
colnames(TP3all)[c(4,13,16)]=c("Plot","GE3","GI3scale")
colnames(GGS_TP1)[c(1,5)]=c("Plot","Env"); colnames(GGS_TP4)[c(1,5)]=c("Plot","Env");colnames(GGS_TP6)[c(1,5)]=c("Plot","Env")
GGS_TP1$Location=paste(GGS_TP1$Location,"20",sep="");GGS_TP4$Location=paste(GGS_TP4$Location,"20",sep="");GGS_TP6$Location=paste(GGS_TP6$Location,"20",sep="");

TP1_1920=rbind(TP1all[,c(2,4,17,12,13,16)], GGS_TP1[,c(4,1,5,7,16,20)])
TP2_1920=rbind(TP2all[,c(2,4,17,12,13,16)], GGS_TP4[,c(4,1,5,7,16,20)])
TP3_1920=rbind(TP3all[,c(2,4,17,12,13,16)], GGS_TP6[,c(4,1,5,7,16,20)])
TP1_1920$Plot=paste(TP1_1920$Plot,"_",TP1_1920$Env, sep=""); TP2_1920$Plot=paste(TP2_1920$Plot,"_",TP2_1920$Env, sep="");TP3_1920$Plot=paste(TP3_1920$Plot,"_",TP3_1920$Env, sep="");

TP1_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif

TP1_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif

#combined BLUE summary for 2019/2020 data
PM6_GI3_1920_BLUE=get_entry_fixed_Effects(TP1_GI3_1920.lm, 'GI3')
PM47_GI3_1920_BLUE=get_entry_fixed_Effects(TP2_GI3_1920.lm, 'GI3')
PM111_GI3_1920_BLUE=get_entry_fixed_Effects(TP3_GI3_1920.lm, 'GI3')
PM6_GE3_1920_BLUE=get_entry_fixed_Effects(TP1_GE3_1920.lm, 'GE3')
PM47_GE3_1920_BLUE=get_entry_fixed_Effects(TP2_GE3_1920.lm, 'GE3')
PM111_GE3_1920_BLUE=get_entry_fixed_Effects(TP3_GE3_1920.lm, 'GE3')

#combined BLUP summary for 2019/2020 data
PM6_GE3_1920_BLUP=gapit_blup(TP1_GE3_1920.lm, 'GE3')
PM47_GE3_1920_BLUP=gapit_blup(TP2_GE3_1920.lm, 'GE3')
PM111_GE3_1920_BLUP=gapit_blup(TP3_GE3_1920.lm, 'GE3')
PM6_GI3_1920_BLUP=gapit_blup(TP1_GI3_1920.lm, 'GI3')
PM47_GI3_1920_BLUP=gapit_blup(TP2_GI3_1920.lm, 'GI3')
PM111_GI3_1920_BLUP=gapit_blup(TP3_GI3_1920.lm, 'GI3')


##### load and process preharvest sprouting data #####

load("PhenotypeData/RawData/allGGSPHSdata.RData")
phs.lm=lmer(PHS ~  Env + (1|Env:Harvest) + Entry +(1|Env:Rep) , data=allGGSphs) #GGS19, GGS20

get_entry_fixed_Effects_phs = function(model, trait){
  x = data.frame(fixef(model))
  x$names = rownames(x)
  x$names2 = substr(x$names,1,5)
  x = subset(x, x$names2 %in% c('(Inte', 'Entry'))
  x$taxa = substring(x$names, 6)
  x$taxa[1] = as.character(sort(unique(model@frame$Entry))[1])
  #intercept = x$coef.model.[1]
  #x$coef.model.[1] = 0
  x = data.frame(taxa = x$taxa, x[,1])
  colnames(x)[2] = trait
  return(x)
} #extract Entry BLUEs

PHS.blues=get_entry_fixed_Effects_phs(phs.lm, "PHS")
PHS.blues=PHS.blues[which(PHS.blues$taxa %in% myGD20_prune$taxa),]


##### run GAPIT models #####

library(GAPIT3)

## use 2 PCs and MAF of 0.05 
##also omit six-rows 
six=c("DH133529", "DH133535", "MS1054111_01", "MS1054115_03", "Tamalpais", "BB5","BB28","BBB528","Purple_Valley","Morex","Steptoe","ND20448","Bohmische_Nackte",
      "Larker","Tamparkorn","Chevron","Bere","Sativum_Jessen_England","Mestny","Odessa","Tiree_6_row","Austrian_Early","Asplund","Kagelkorn","Hatif_de_Grignon",
      "Manchuria","Monte_Cristo","Riojana","Oderbrucker","PI_383933")
#PHS
PHS.MLMM.BLUE.GGS20 =GAPIT(Y= PHS.blues, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2,file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP1 2020
GE3tp1.MLMM.BLUE.GGS20 =GAPIT(Y= TP1_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GE5tp1.MLMM.BLUE.GGS20 =GAPIT(Y= TP1_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp1.MLMM.BLUE.GGS20 =GAPIT(Y= TP1_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5tp1.MLMM.BLUE.GGS20 =GAPIT(Y= TP1_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP2 2020
GE3TP2.MLMM.BLUE.GGS20 =GAPIT(Y= TP2_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GE5TP2.MLMM.BLUE.GGS20 =GAPIT(Y= TP2_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3TP2.MLMM.BLUE.GGS20 =GAPIT(Y= TP2_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5TP2.MLMM.BLUE.GGS20 =GAPIT(Y= TP2_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP3 2020
GE3TP3.MLMM.BLUE.GGS20 =GAPIT(Y= TP3_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GE5TP3.MLMM.BLUE.GGS20 =GAPIT(Y= TP3_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3TP3.MLMM.BLUE.GGS20 =GAPIT(Y= TP3_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5TP3.MLMM.BLUE.GGS20 =GAPIT(Y= TP3_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP4 2020
GE3TP4.MLMM.BLUE.GGS20 =GAPIT(Y= TP4_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GE5TP4.MLMM.BLUE.GGS20 =GAPIT(Y= TP4_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3TP4.MLMM.BLUE.GGS20 =GAPIT(Y= TP4_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5TP4.MLMM.BLUE.GGS20 =GAPIT(Y= TP4_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP5 2020
GE3TP5.MLMM.BLUE.GGS20 =GAPIT(Y= TP5_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GE5TP5.MLMM.BLUE.GGS20 =GAPIT(Y= TP5_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3TP5.MLMM.BLUE.GGS20 =GAPIT(Y= TP5_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5TP5.MLMM.BLUE.GGS20 =GAPIT(Y= TP5_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP6 2020
  #no variance for GE at TP6
#GE3TP6.MLMM.BLUE.GGS20 =GAPIT(Y= TP6_GGS20_BLUE_all[,c(1,2)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
#GE5TP6.MLMM.BLUE.GGS20 =GAPIT(Y= TP6_GGS20_BLUE_all[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3TP6.MLMM.BLUE.GGS20 =GAPIT(Y= TP6_GGS20_BLUE_all[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI5TP6.MLMM.BLUE.GGS20 =GAPIT(Y= TP6_GGS20_BLUE_all[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP1 2019
GE3tp1.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,3)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp1.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,6)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP2 2019
GE3tp2.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,4)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp2.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,7)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP3 2019
GE3tp3.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,5)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp3.MLMM.BLUE.GGS19 =GAPIT(Y= GGS_BLUEs_2019[,c(1,8)], GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP1/PM6 2019+2020
GE3tp1.MLMM.BLUE.GGS1920 =GAPIT(Y= PM6_GE3_1920_BLUE, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp1.MLMM.BLUE.GGS1920 =GAPIT(Y= PM6_GI3_1920_BLUE, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP2/PM47 2019+2020
GE3tp2.MLMM.BLUE.GGS1920 =GAPIT(Y= PM47_GE3_1920_BLUE, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
GI3tp2.MLMM.BLUE.GGS1920 =GAPIT(Y= PM47_GI3_1920_BLUE, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)

#TP3/PM111 2019+2020
GI3tp3.MLMM.BLUE.GGS1920 =GAPIT(Y= PM111_GI3_1920_BLUE, GD=myGD20_prune[-c(which(rownames(myGD20_prune) %in% six)),], GM=myGM20_prune, PCA.total=2, file.output=F, Geno.View.output=F, model="MLMM", Major.allele.zero = F, SNP.MAF=0.05)
