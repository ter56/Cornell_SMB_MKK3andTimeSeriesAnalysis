#GGS20 raw data process
#run the script 'GGS_2020_rawdata_process.R' first.
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


###load 2019 raw pheno data
#IMPORTANT: 2019 TP1 matches with 2020 TP1, 2019 TP2 matches with 2020 TP4, and 2019 TP3 matches with 2020 TP6
load("PhenotypeData/RawData/2019/TP1 all data for analysis.RData")
load("PhenotypeData/RawData/2019/TP2 all data for analysis.RData")
load("PhenotypeData/RawData/2019/TP3 30k all data for analysis.RData");TP3all=TP3all_30k; rm(TP3all_30k)

#address outliers
GGS_TP1[152,16:21] = NA #incomplete data entry
GGS_TP5[c(2204,904),16:21] = NA #very low germ, likely due to disease
GGS_TP6[which(GGS_TP6$KernelsLeft>5),16:21] =NA #very low germ, likely due to disease
GGS_TP7[which(GGS_TP7$KernelsLeft>5),16:21] =NA #very low germ, likely due to disease and poor seed quality at the bottom of the bag...

hist(GGS_TP7$GE3)
#heritability estimates 
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
  #TP7
  Cullis_H2(lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),])) # .665
  Cullis_H2(lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),])) #.692
  Cullis_H2(lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),])) #.855
  Cullis_H2(lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),])) #.855
  
#fixed effects models
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

tp7_ge5.lm=lm(GE5~Location  + rep  + Entry, data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
tp7_ge3.lm=lm(GE3~Location  + rep  + Entry, data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
tp7_gi5.lm=lm(GI5scale~Location  + rep  + Entry, data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif
tp7_gi3.lm=lm(GI3scale~Location  + rep  + Entry, data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif

#random effects models
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

tp7_ge5.lmer=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
tp7_ge3.lmer=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
tp7_gi5.lmer=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif
tp7_gi3.lmer=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif

## BLUEs for each trait, organized by time point
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
}
get_entry_blups=function(mod,trait){
  df=ranef(mod)$'Entry' + fixef(mod)[1]
  df=cbind(rownames(df), df); 
  colnames(df)=c("taxa", trait); 
  df[,1]=as.character(df[,1])
  df
}

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
}
BLUP_summary=function(mods,PM,TP,date){
  l=list(get_entry_blups(mods[[1]], 'GE3'),get_entry_blups(mods[[2]], 'GE5'),
         get_entry_blups(mods[[3]], 'GI3scale'),get_entry_blups(mods[[4]], 'GI5scale'))
  
  TP_summary <- Reduce(
    function(x, y, ...) merge(x, y, by="taxa" ,all = TRUE),
    l
  )
  TP_summary=TP_summary[-c(which(!TP_summary$taxa %in% myGD20_prune$taxa)),]
  TP_summary$PM_date=PM
  TP_summary$TP=TP
  TP_summary$date=date
  TP_summary
}


setwd('PhenotypeData/ProcessedData/2020/')

TP1_GGS20_BLUE_all=BLUE_summary(list(tp1_ge3.lm,tp1_ge5.lm,tp1_gi3.lm,tp1_gi5.lm), PM=5, TP="TP1", date="9/27/2020")
TP2_GGS20_BLUE_all=BLUE_summary(list(tp2_ge3.lm,tp2_ge5.lm,tp2_gi3.lm,tp2_gi5.lm), PM=19, TP="TP2", date="10/19/2020")
TP3_GGS20_BLUE_all=BLUE_summary(list(tp3_ge3.lm,tp3_ge5.lm,tp3_gi3.lm,tp3_gi5.lm), PM=33, TP="TP3", date="10/26/2020")
TP4_GGS20_BLUE_all=BLUE_summary(list(tp4_ge3.lm,tp4_ge5.lm,tp4_gi3.lm,tp4_gi5.lm), PM=47, TP="TP4", date="11/9/2020")
TP5_GGS20_BLUE_all=BLUE_summary(list(tp5_ge3.lm,tp5_ge5.lm,tp5_gi3.lm,tp5_gi5.lm), PM=68, TP="TP5", date="11/30/2020")
TP6_GGS20_BLUE_all=BLUE_summary(list(tp6_ge3.lm,tp6_ge5.lm,tp6_gi3.lm,tp6_gi5.lm), PM=110, TP="TP6", date="1/11/2021")
TP7_GGS20_BLUE_all=BLUE_summary(list(tp7_ge3.lm,tp7_ge5.lm,tp7_gi3.lm,tp7_gi5.lm), PM=159, TP="TP7", date="3/1/2021")

all_BLUE = rbind(TP1_GGS20_BLUE_all,
                 TP2_GGS20_BLUE_all,
                 TP3_GGS20_BLUE_all,
                 TP4_GGS20_BLUE_all,
                 TP5_GGS20_BLUE_all,
                 TP6_GGS20_BLUE_all,
                 TP7_GGS20_BLUE_all)


save(TP1_GGS20_BLUE_all, file="TP1_GGS20_BLUE_all.RData")
save(TP2_GGS20_BLUE_all, file="TP2_GGS20_BLUE_all.RData")
save(TP3_GGS20_BLUE_all, file="TP3_GGS20_BLUE_all.RData")
save(TP4_GGS20_BLUE_all, file="TP4_GGS20_BLUE_all.RData")
save(TP5_GGS20_BLUE_all, file="TP5_GGS20_BLUE_all.RData")
save(TP6_GGS20_BLUE_all, file="TP6_GGS20_BLUE_all.RData")
save(TP7_GGS20_BLUE_all, file="TP7_GGS20_BLUE_all.RData")
save(all_BLUE, file = "GGS2020_BLUE_summary_allTP.RData")

TP1_GGS20_BLUP_all=BLUP_summary(list(tp1_ge3.lmer,tp1_ge5.lmer,tp1_gi3.lmer,tp1_gi5.lmer), PM=5, TP="TP1", date="9/27/2020")
TP2_GGS20_BLUP_all=BLUP_summary(list(tp2_ge3.lmer,tp2_ge5.lmer,tp2_gi3.lmer,tp2_gi5.lmer), PM=19, TP="TP2", date="10/19/2020")
TP3_GGS20_BLUP_all=BLUP_summary(list(tp3_ge3.lmer,tp3_ge5.lmer,tp3_gi3.lmer,tp3_gi5.lmer), PM=33, TP="TP3", date="10/26/2020")
TP4_GGS20_BLUP_all=BLUP_summary(list(tp4_ge3.lmer,tp4_ge5.lmer,tp4_gi3.lmer,tp4_gi5.lmer), PM=47, TP="TP4", date="11/9/2020")
TP5_GGS20_BLUP_all=BLUP_summary(list(tp5_ge3.lmer,tp5_ge5.lmer,tp5_gi3.lmer,tp5_gi5.lmer), PM=68, TP="TP5", date="11/30/2020")
TP6_GGS20_BLUP_all=BLUP_summary(list(tp6_ge3.lmer,tp6_ge5.lmer,tp6_gi3.lmer,tp6_gi5.lmer), PM=110, TP="TP6", date="1/11/2021")
TP7_GGS20_BLUP_all=BLUP_summary(list(tp7_ge3.lmer,tp7_ge5.lmer,tp7_gi3.lmer,tp7_gi5.lmer), PM=159, TP="TP7", date="3/1/2021")

all_BLUP = rbind(TP1_GGS20_BLUP_all,
                 TP2_GGS20_BLUP_all,
                 TP3_GGS20_BLUP_all,
                 TP4_GGS20_BLUP_all,
                 TP5_GGS20_BLUP_all,
                 TP6_GGS20_BLUP_all,
                 TP7_GGS20_BLUP_all)

save(TP1_GGS20_BLUP_all, file="TP1_GGS20_BLUP_all.RData")
save(TP2_GGS20_BLUP_all, file="TP2_GGS20_BLUP_all.RData")
save(TP3_GGS20_BLUP_all, file="TP3_GGS20_BLUP_all.RData")
save(TP4_GGS20_BLUP_all, file="TP4_GGS20_BLUP_all.RData")
save(TP5_GGS20_BLUP_all, file="TP5_GGS20_BLUP_all.RData")
save(TP6_GGS20_BLUP_all, file="TP6_GGS20_BLUP_all.RData")
save(TP7_GGS20_BLUP_all, file="TP7_GGS20_BLUP_all.RData")

save(all_BLUP, file = 'GGS2020_BLUP_summary_allTP.Rdata')

setwd(rprojroot::find_rstudio_root_file())
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

TP1all$Location=ifelse(substr(TP1all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP2all$Location=ifelse(substr(TP2all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP3all$Location=ifelse(substr(TP3all$PLOT,1,1)=="1", "Caldwell19","Ketola19")
TP1all$rep=ifelse(TP1all$rep=="A", "1","2")
TP2all$rep=ifelse(TP2all$rep=="A", "1","2")
TP3all$rep=ifelse(TP3all$rep=="A", "1","2")

colnames(TP1all)[c(13,16)]=c("GE3","GI3scale")
colnames(TP2all)[c(13,16)]=c("GE3","GI3scale")
colnames(TP3all)[c(13,16)]=c("GE3","GI3scale")

TP1_1920=rbind(TP1all[,c(2,17,12,13,16)], GGS_TP1[,c(4,5,7,16,20)])
TP2_1920=rbind(TP2all[,c(2,17,12,13,16)], GGS_TP4[,c(4,5,7,16,20)])
TP3_1920=rbind(TP3all[,c(2,17,12,13,16)], GGS_TP6[,c(4,5,7,16,20)])


TP1_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GE3_1920.lm=lm(GE3~Location  + Location:rep  + Entry, data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif

TP1_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GI3_1920.lm=lm(GI3scale~Location  + Location:rep  + Entry, data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif

TP1_GE3_1920.lmer=lmer(GE3~Location  + Location:rep  + (1|Entry), data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GE3_1920.lmer=lmer(GE3~Location  + Location:rep  + (1|Entry), data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GE3_1920.lmer=lmer(GE3~Location  + Location:rep  + (1|Entry), data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif

TP1_GI3_1920.lmer=lmer(GI3scale~Location  + Location:rep  + (1|Entry), data=TP1_1920[-c(which(TP1_1920$Location=="check")),]) #all fixef are signif
TP2_GI3_1920.lmer=lmer(GI3scale~Location  + Location:rep  + (1|Entry), data=TP2_1920[-c(which(TP2_1920$Location=="check")),]) #all fixef are signif
TP3_GI3_1920.lmer=lmer(GI3scale~Location  + Location:rep  + (1|Entry), data=TP3_1920[-c(which(TP3_1920$Location=="check")),]) #all fixef are signif


#combined BLUP summary for 2019/2020 data
TP1_GE3_1920_BLUP=get_entry_blups(TP1_GE3_1920.lmer, 'GE3')
TP2_GE3_1920_BLUP=get_entry_blups(TP2_GE3_1920.lmer, 'GE3')
TP3_GE3_1920_BLUP=get_entry_blups(TP3_GE3_1920.lmer, 'GE3')
TP1_GI3_1920_BLUP=get_entry_blups(TP1_GI3_1920.lmer, 'GI3')
TP2_GI3_1920_BLUP=get_entry_blups(TP2_GI3_1920.lmer, 'GI3')
TP3_GI3_1920_BLUP=get_entry_blups(TP3_GI3_1920.lmer, 'GI3')

#combined BLUE summary for 2019/2020 data
PM6_GI3_1920_BLUE=get_entry_fixed_Effects(TP1_GI3_1920.lm, 'GI3')
PM47_GI3_1920_BLUE=get_entry_fixed_Effects(TP2_GI3_1920.lm, 'GI3')
PM111_GI3_1920_BLUE=get_entry_fixed_Effects(TP3_GI3_1920.lm, 'GI3')
PM6_GE3_1920_BLUE=get_entry_fixed_Effects(TP1_GE3_1920.lm, 'GE3')
PM47_GE3_1920_BLUE=get_entry_fixed_Effects(TP2_GE3_1920.lm, 'GE3')
PM111_GE3_1920_BLUE=get_entry_fixed_Effects(TP3_GE3_1920.lm, 'GE3')

setwd('PhenotypeData/ProcessedData/2019and2020Combined/')
# save(PM6_GI3_1920_BLUE, file="PM6_GI3_1920_BLUE.RData")
# save(PM47_GI3_1920_BLUE, file="PM47_GI3_1920_BLUE.RData")
# save(PM111_GI3_1920_BLUE, file="PM111_GI3_1920_BLUE.RData")
# save(PM6_GE3_1920_BLUE, file="PM6_GE3_1920_BLUE.RData")
# save(PM47_GE3_1920_BLUE, file="PM47_GE3_1920_BLUE.RData")
# save(PM111_GE3_1920_BLUE, file="PM111_GE3_1920_BLUE.RData")


##### load and process preharvest sprouting data #####
setwd(rprojroot::find_rstudio_root_file())
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



#summaries
get_entry_blups=function(mod,trait){
  df=ranef(mod)$'Entry' + fixef(mod)[1]
  df=cbind(rownames(df), df); 
  colnames(df)=c("taxa", trait); 
  df[,1]=as.character(df[,1])
  df
}

q=list(get_entry_blups(tp1_ge3.lmer, "GE3"),get_entry_blups(tp2_ge3.lmer, "GE3"),get_entry_blups(tp3_ge3.lmer, "GE3"),
       get_entry_blups(tp4_ge3.lmer, "GE3"),get_entry_blups(tp5_ge3.lmer, "GE3"),get_entry_blups(tp6_ge3.lmer, "GE3") )

ge3_summary <- Reduce(
  function(x, y, ...) merge(x, y, by="taxa" ,all = TRUE),
  q
)
ge3_summary=ge3_summary[-c(which(!ge3_summary$taxa %in% myGD20$taxa)),]
colnames(ge3_summary)[2:6] = c("ge3_TP1","ge3_TP2","ge3_TP3","ge3_TP4","ge3_TP5")

ge3_summary$AlaAT=round(myGD20[which(myGD20$taxa %in% ge3_summary$taxa),which(colnames(myGD20)=="Qsd1")],0)
ge3_summary$GA20ox1=round(myGD20[which(myGD20$taxa %in% ge3_summary$taxa),which(colnames(myGD20)=="SCRI_RS_99344")],0)
ge3_summary$MKK3=round(myGD20[which(myGD20$taxa %in% ge3_summary$taxa),which(colnames(myGD20)=="JHI-Hv50k-2016-367342")],0)
ge3_summary$haplo=as.factor(paste(ge3_summary$AlaAT, ge3_summary$GA20ox1, ge3_summary$MKK3, sep=""))

GGSmelt=melt(ge3_summary, id.vars =c("taxa","haplo"), measure.vars = c( "ge3_TP1", "ge3_TP2","ge3_TP3", "ge3_TP4","ge3_TP5"))
GGSmelt=GGSmelt[-c(which(GGSmelt$haplo %in% c("012","021","102","120","122","212"))),]
GGSmelt=droplevels(GGSmelt)
GGSmelt$haplo=factor(GGSmelt$haplo, levels = c("202", "002","220", "020","222", "022" ,"200","000"))

ggplot(data=droplevels(GGSmelt), aes(x=haplo, y=value, fill=haplo)) +
  geom_boxplot()+
  theme_bw() + ylim(0,1)+
  facet_wrap(~variable, nrow=1) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "202") +
  scale_x_discrete(breaks=c("202","002", "220","020", "222","022", "000","200"), labels=c("DNN","NNN", "DDD","NDD","DDN", "NDN", "NND","DND")) +
  #ggtitle("CU1 PHS")+ 
  theme(legend.position = "none")+
 # scale_fill_manual(values=c("#b2182b","#ef8a62","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac")) +
  xlab("") + ylab("Germination energy") 
