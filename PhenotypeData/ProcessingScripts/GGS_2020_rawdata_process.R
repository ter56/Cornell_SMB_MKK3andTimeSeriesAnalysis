# input 2020 GGS germination data
# create GE and GI phenotypes, check distribution, look at check lines and coolers
# consolidate files into one? At least make R objects to read in faster. 

library(here)
library(readxl)
library(lme4)

setwd(rprojroot::find_rstudio_root_file())

GGS_TP1=read_excel("PhenotypeData/RawData/2020/GGSTp1_all.xlsx")
GGS_TP2=read_excel("PhenotypeData/RawData/2020/GGSTp2_all.xlsx")
GGS_TP3=read_excel("PhenotypeData/RawData/2020/GGSTp3_all.xlsx")
GGS_TP4=read_excel("PhenotypeData/RawData/2020/GGSTp4_all.xlsx")
GGS_TP5=read_excel("PhenotypeData/RawData/2020/GGSTp5_all.xlsx")
GGS_TP6=read_excel("PhenotypeData/RawData/2020/GGSTp6_all.xlsx")
GGS_TP7=read_excel("PhenotypeData/RawData/2020/GGSTp7_all.xlsx")

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
}


GGS_TP1=germ_phenos(GGS_TP1)
GGS_TP2=germ_phenos(GGS_TP2)
GGS_TP3=germ_phenos(GGS_TP3)
GGS_TP4=germ_phenos(GGS_TP4)
GGS_TP5=germ_phenos(GGS_TP5)
GGS_TP6=germ_phenos(GGS_TP6)
GGS_TP7=germ_phenos(GGS_TP7)


setwd('PhenotypeData/RawData/2020/')
# save(GGS_TP1, file="GGS_TP1_for_analysis.RData")
# save(GGS_TP2, file="GGS_TP2_for_analysis.RData")
# save(GGS_TP3, file="GGS_TP3_for_analysis.RData")
# save(GGS_TP4, file="GGS_TP4_for_analysis.RData")
# save(GGS_TP5, file="GGS_TP5_for_analysis.RData")
# save(GGS_TP6, file="GGS_TP6_for_analysis.RData")
# save(GGS_TP7, file="GGS_TP7_for_analysis.RData")

setwd(rprojroot::find_rstudio_root_file())


#are cooler effects significant?
tp1_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif
tp1_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #all fixef are signif
tp1_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif
tp1_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP1[-c(which(GGS_TP1$Location=="check")),]) #rep not signif

tp2_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
#tp2_ge5.lm2=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #all fixef are signif
tp2_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif
tp2_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP2[-c(which(GGS_TP2$Location=="check")),]) #rep not signif

tp3_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
#tp3_gi5.lm2=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif
tp3_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP3[-c(which(GGS_TP3$Location=="check")),]) #all fixef are signif

tp4_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
#tp4_gi5.lm2=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif
tp4_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP4[-c(which(GGS_TP4$Location=="check")),]) #all fixef are signif

tp5_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #all fixef are signif
tp5_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif
tp5_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP5[-c(which(GGS_TP5$Location=="check")),]) #rep not signif

tp6_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #all fixef are signif
tp6_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif
tp6_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP6[-c(which(GGS_TP6$Location=="check")),]) #rep not signif

tp7_ge5.lm=lmer(GE5~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
tp7_ge3.lm=lmer(GE3~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #Rep not signif
tp7_gi5.lm=lmer(GI5scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif
tp7_gi3.lm=lmer(GI3scale~Location  + rep  + (1|Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #rep not signif

tp7_ge5.lm1=lm(GE5~Location  + rep  + (Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
anova(tp7_ge5.lm1)
tp7_ge3.lm1=lm(GE3~Location  + rep  + (Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef are signif
anova(tp7_ge3.lm1)
tp7_gi5.lm1=lm(GI5scale~Location  + rep  + (Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef signif, location marginally
anova(tp7_gi5.lm1)
tp7_gi3.lm1=lm(GI3scale~Location  + rep  + (Entry), data=GGS_TP7[-c(which(GGS_TP7$Location=="check")),]) #all fixef signif, location only marginally
anova(tp7_gi3.lm1)

#residual plots look fine for all models 


GGS20 = data.frame(read_excel("PhenotypeData/RawData/2020/GGS_20 fieldbook_DS.xlsx"))
GGS20$Entry[440]="Megs_song"
GGS20$Entry=gsub("\\.","_",GGS20$Entry)
GGS20$Entry=gsub("\\-","_",GGS20$Entry)
GGS20$Entry=gsub(" ","_",GGS20$Entry)


blups_for_hist = function(mod, TP){
  df=as.data.frame(ranef(mod)$Entry)
  colnames(df)[1]="blup"
  df$blup=df$blup+fixef(mod)[1]
  df$Entry=rownames(df)
  df=merge(df, GGS20[,1:3], by="Entry", all=TRUE)
  df=df[-c(which(df$Population %in% "check")),]
  df$TP=as.factor(TP)
  df
}

GE3=rbind(blups_for_hist(tp1_ge3.lm, "TP1"), blups_for_hist(tp2_ge3.lm, "TP2"),blups_for_hist(tp3_ge3.lm, "TP3"),
          blups_for_hist(tp4_ge3.lm, "TP4"),blups_for_hist(tp5_ge3.lm, "TP5"),blups_for_hist(tp6_ge3.lm, "TP6"),
          blups_for_hist(tp7_ge3.lm, "TP7"))

GE5=rbind(blups_for_hist(tp1_ge5.lm, "TP1"), blups_for_hist(tp2_ge5.lm, "TP2"),blups_for_hist(tp3_ge5.lm, "TP3"),
          blups_for_hist(tp4_ge5.lm, "TP4"),blups_for_hist(tp5_ge5.lm, "TP5"),blups_for_hist(tp6_ge5.lm, "TP6"),
          blups_for_hist(tp7_ge5.lm, 'TP7'))

GI3=rbind(blups_for_hist(tp1_gi3.lm, "TP1"), blups_for_hist(tp2_gi3.lm, "TP2"),blups_for_hist(tp3_gi3.lm, "TP3"),
          blups_for_hist(tp4_gi3.lm, "TP4"),blups_for_hist(tp5_gi3.lm, "TP5"),blups_for_hist(tp6_gi3.lm, "TP6"),
          blups_for_hist(tp7_gi3.lm, 'TP7'))

GI5=rbind(blups_for_hist(tp1_gi5.lm, "TP1"), blups_for_hist(tp2_gi5.lm, "TP2"),blups_for_hist(tp3_gi5.lm, "TP3"),
          blups_for_hist(tp4_gi5.lm, "TP4"),blups_for_hist(tp5_gi5.lm, "TP5"),blups_for_hist(tp6_gi5.lm, "TP6"),
          blups_for_hist(tp7_gi5.lm, 'TP7'))


blups_for_hist2 = function(mod, TP){
  df=as.data.frame(ranef(mod)$Entry)
  colnames(df)[1]="blup"
  df$blup=df$blup+fixef(mod)[1]
  df$Entry=rownames(df)
  df=merge(df, GGS20[,1:3], by="Entry", all=TRUE)
  df$TP=as.factor(TP)
  df
}
GI3row=data.frame(blups_for_hist2(tp1_gi3.lm, "TP1")[,c(1,3,4)],
          TP1 = blups_for_hist2(tp1_gi3.lm, "TP1")[,2], 
          TP2 = blups_for_hist2(tp2_gi3.lm, "TP2")[,2],
          TP3 = blups_for_hist2(tp3_gi3.lm, "TP3")[,2],
          TP4 = blups_for_hist2(tp4_gi3.lm, "TP4")[,2],
          TP5 = blups_for_hist2(tp5_gi3.lm, "TP5")[,2],
          TP6 = blups_for_hist2(tp6_gi3.lm, "TP6")[,2])
GI3withchecks = rbind(blups_for_hist2(tp1_gi3.lm, "TP1"), blups_for_hist2(tp2_gi3.lm, "TP2"),blups_for_hist2(tp3_gi3.lm, "TP3"),blups_for_hist2(tp4_gi3.lm, "TP4"),blups_for_hist2(tp5_gi3.lm, "TP5"),blups_for_hist2(tp6_gi3.lm, "TP6"))


GE3$blup[which(GE3$blup>1)] =1
GE5$blup[which(GE5$blup>1)] =1

library(ggplot2)

ggplot(GE3) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 3 day GE distribution")

ggplot(GE5) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 5 day GE distribution")

ggplot(GI3) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 3 day GI distribution")

ggplot(GI5) +
  geom_density( aes(x=blup, group=TP, fill=TP),alpha=0.5) +
  facet_wrap(~Population, nrow=5) +
  theme_bw()+
  ggtitle("GGS 5 day GI distribution")

aggregate(blup ~ TP +Population, GI3, mean)
aggregate(blup ~ TP +Population, GI5, mean)


#checks
GGSall=rbind(GGS_TP1, GGS_TP2[,-c(16)], GGS_TP3, GGS_TP4, GGS_TP5, GGS_TP6, GGS_TP7)
GGSall$TP=rep(1:7, each=2208)
aggregate(GE3 ~ TP + Entry, GGSall[which(GGSall$Location=="check"),], mean)
aggregate(GE5 ~ TP + Entry, GGSall[which(GGSall$Location=="check"),], mean)

aggregate(GI3 ~ TP + Entry, GGSall[which(GGSall$Location=="check"),], mean)
ggplot(aggregate(GI3 ~ TP + Entry, GGSall[which(GGSall$Location=="check"),], mean),aes(TP,GI3,color =Entry))+geom_point()
aggregate(GI5 ~ TP + Entry, GGSall[which(GGSall$Location=="check"),], mean)

checks = GGSall %>% filter(Location == 'check') %>% group_by(TP) %>% summarise(mean = mean(GI3)) %>% mutate(meandev = mean -mean[1],
                                                                                                            time = c(0,14,28,42,63,105,154))

#there do seem to be differences between TP for checks, especially for GI Unless that freezer does not keep them cold enough...
# Sum summary plots for a small presentation ####

ggplot(GE3) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 3 day GE distribution")

ggplot(GE5) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 5 day GE distribution")

ggplot(GI3) +
  geom_density( aes(x=blup, group=Population, fill=Population),alpha=0.5) +
  facet_wrap(~TP) +
  theme_bw()+
  ggtitle("GGS 3 day GI distribution")

ggplot(GI5) +
  geom_density( aes(x=blup, group=TP, fill=TP),alpha=0.5) +
  facet_wrap(~Population, nrow=5) +
  theme_bw()+
  ggtitle("GGS 5 day GI distribution")

ggplot(GI5, aes(x = TP, y = blup, group=Entry))+geom_line()




