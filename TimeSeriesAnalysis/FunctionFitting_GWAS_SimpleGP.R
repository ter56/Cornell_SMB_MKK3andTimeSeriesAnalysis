library(Hmisc);library(drc);library(ggplot2);library(readxl);library(reshape2);
library(patchwork);library(rrBLUP);library(plyr);library(dplyr);
library(knitr);library(tidyr); library(fda) ; library(magic) 
library(patchwork); library(knitr);library(stringr);library(LDcorSV)
library(ggtext);library(ggh4x)
library(GAPIT3)

setwd(rprojroot::find_rstudio_root_file())
source('TimeSeriesAnalysis/Functions/FPCA_function.R')
source('TimeSeriesAnalysis/Functions/pca_fun.R')
source('TimeSeriesAnalysis/Functions/pca_score.R')
source('TimeSeriesAnalysis/Functions/tuning_nointer.R')
DF_OperationsV2 = function(dataframe){
  dataframe = dataframe[order(dataframe$P.value),]
  dataframe$logPercentileQQplot = -log10(c(1:length(dataframe$SNP))/length(dataframe$SNP))
  dataframe$rank = c(1:length(dataframe$SNP))
  dataframe$FDRPval = length(dataframe$SNP)*dataframe$P.value/dataframe$rank
  dataframe = dataframe[order(dataframe$Chromosome, as.numeric(dataframe$Position)),]
  dataframe$log10PVal = -log10(dataframe$P.value)
  dataframe$ordinal = c(1:length(dataframe$SNP))
  return(dataframe)
} #Operations to be performed on GAPIToutput$GWAS so that homeMadeManhattanLDPruned works
TableOutput = function(GWAS.sum.dataframe, n = 20){
  GWAS.sum.dataframe[order(GWAS.sum.dataframe$rank),c('SNP','Chromosome','Position','P.value','maf','FDRPval')][1:n,] %>% 
    mutate(maf = as.numeric(as.character(maf))) %>%
    kable(digits = c(0,1,0,8,2,6), align = 'lcccrr')}
ld_heatmap=function(df, markerList){
  ld <- as.matrix(round(df,0))
  
  if(c(-1,3,4) %in% ld){
    ld[which(ld==3)]=2
    ld[which(ld==4)]=2
    ld[which(ld== -1)]=0
  }
  
  LD <- LD.Measures(donnees=ld,  na.presence=F)
  #LD$loc1=as.character(LD$loc1); LD$loc2=as.character(LD$loc2)
  r2 <- matrix(0, nrow=ncol(df), ncol=ncol(df))
  r2[lower.tri(r2, diag=FALSE)] <- LD$r2
  r2 <- t(r2)
  r2 <- as.data.frame(round(r2, 5))
  diag(r2) <- 1
  r2[lower.tri(r2)] = NA
  rownames(r2)=colnames(df); colnames(r2)=rownames(r2)
  r_2=melt(as.matrix(r2), na.rm=T)
  
  r_2  = r_2 %>% 
    mutate(ChrPos = mapvalues(as.character(Var1), from = as.character(markerList$SNP),
                              to = paste0(markerList$Chromosome,':',
                                          str_pad(as.character(markerList$Position),pad = "0",width = 9,side = 'left'))))
  
  graphic = ggplot(r_2, aes(Var2, ChrPos, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", name="r2") +
    theme_classic() #+    ggtitle(paste("LD r2 from",colnames(r2)[1],"-", colnames(r2)[length(colnames(r2))], sep=" " ))
  return(graphic)
}

# Loading the BLUES from all lines and GAPIT functions #####
setwd(rprojroot::find_rstudio_root_file())
load('PhenotypeData/ProcessedData/2020/GGS2020_BLUE_summary_allTP.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM111_GE3_1920_BLUE.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM47_GE3_1920_BLUE.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM6_GE3_1920_BLUE.RData')
load('GenotypeData/myGD_LDpruned_w_KASP.RData')
load('GenotypeData/myGM_LDpruned_w_KASP.RData')

# Bind BLUEs Together and change results that are above 1 to 1 as noted in script below#####
All.BluesGE31920 = all_BLUE  %>% mutate(time = PM_date-5) %>% filter(TP %in% c('TP2','TP3','TP5', 'TP7')) %>% 
  select(taxa, GE3, TP, time) %>% 
  rbind(.,PM111_GE3_1920_BLUE %>% mutate(TP = 'TP6', time = 105),
        PM47_GE3_1920_BLUE %>% mutate(TP = 'TP4', time = 42),
        PM6_GE3_1920_BLUE %>% mutate(TP = 'TP1', time = 0)) %>% arrange(taxa, TP) %>%
  mutate(GE3 = ifelse(GE3>0,ifelse(GE3<1,GE3,1),0))

All.BluesGE31920 %>% group_by(taxa) %>% dplyr::tally() %>% select(n) %>%table() 
GE3GT4Obs = All.BluesGE31920 %>% group_by(taxa) %>% dplyr::tally() %>% filter(n>4) %>% select(taxa) %>% 
  mutate(taxa = unique(taxa)) %>% ungroup()

# Fit normal logistic models with 3 parameters for all the lines that we can ####
# Three parameter logistic function is used as upper bound is fixed - ie at 100% germ:
GE3_3P_logFits = All.BluesGE31920 %>% 
  filter(taxa %in% GE3GT4Obs$taxa & 
           taxa %nin% c('Chevallier_D10',"G_31_5","P_16_1",'P_29_5')) %>%
  arrange(taxa, TP) %>% group_by(taxa) %>% 
  group_modify(~ broom::tidy(drm(GE3~time, data = .x,
                                 fct=LL.4(fixed = c(NA, NA, 1, NA),
                                          names = c('Rate','Lower','Upper','Centering'))))) %>%
  do(add_row(., taxa = .$taxa[1],term = 'TimeTo90',curve ='Derived',
             estimate = as.double(exp(log((1-.[2,4])/(0.90-.[2,4])-1)/.[1,4]+log(.[3,4]))))%>%
       add_row(., taxa = .$taxa[1],term = 'TimeTo95',curve ='Derived',
               estimate = as.double(exp(log((1-.[2,4])/(0.95-.[2,4])-1)/.[1,4]+log(.[3,4])))) %>%
       add_row(., taxa = .$taxa[1],term = 'rTimeTo95',curve ='Derived',
               estimate = as.double(.[3,4]*exp(log((1-0.95)/0.95)/.[1,4]))) %>%
       add_row(., taxa = .$taxa[1],term = 'rTimeTo90',curve ='Derived',
               estimate = as.double(.[3,4]*exp(log((1-0.90)/0.90)/.[1,4])))
  ) %>%
  ungroup() %>% mutate(estimate = as.numeric(ifelse(is.nan(estimate),0, estimate))) 
GE3_3P_logFits %>% ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')

GE3logTaxatoFilter = GE3_3P_logFits %>% filter(term == 'TimeTo90' & estimate >165 |
                                       term == 'Centering' & estimate>150 | 
                                       (term == 'Rate' & (estimate > 0 | estimate < -13)))
GE3_3P_logFits %>% group_by(term) %>% summarize(mean = mean(estimate), stand.dev = sd(estimate), median = median(estimate))
GE3_3P_logFits %>% filter(taxa %nin% GE3logTaxatoFilter$taxa) %>%group_by(term) %>% summarize(mean = mean(estimate), 
                                                                                    stand.dev = sd(estimate), median = median(estimate),
                                                                                    max = max(estimate),
                                                                                    min = min(estimate))
unique(GE3logTaxatoFilter$taxa) #filters these lines, as well as the problematic lines above: 27 lines removed. 
c('Chevallier_D10',"G_31_5","P_16_1",'P_29_5')
GE3_3P_logFits %>% filter(taxa %nin% GE3logTaxatoFilter$taxa) %>%
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')


GE3fittedCurves = data.frame()
for (i in unique(GE3_3P_logFits$taxa)){
  time = seq(0,150,3)
  tmp = GE3_3P_logFits %>% filter(taxa ==i)
  y = tmp$estimate[2]+(1-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[3]))))
  GE3fittedCurves = rbind(GE3fittedCurves, data.frame(taxa = i, time = time, GE3estimate = y))
}
FacetsParams = GE3_3P_logFits %>% filter(taxa %nin% GE3logTaxatoFilter$taxa) %>%
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')

GE3_3P_logFits %>% filter(taxa %nin% GE3logTaxatoFilter$taxa) %>% group_by(taxa) %>% #can we color based on the time?
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')

#Now we have our variables - lets run the MLMM models on them and see if we need to filter things out.
GE3_3P_logFits.W = GE3_3P_logFits %>% filter(taxa %nin% GE3logTaxatoFilter$taxa) %>% 
  select(term, taxa, estimate) %>% pivot_wider(names_from = 'term', values_from = 'estimate')%>% 
  mutate(rTimeTo90 = ifelse(rTimeTo90>200,NA,rTimeTo90),
         rTimeTo95 = ifelse(rTimeTo95>250,NA,rTimeTo95))
  
dim(GE3_3P_logFits.W) #n = 465

panel.cor <- function(x, y){usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- round(cor(x, y, use = 'complete.obs'), digits=2)
txt <- paste0("R = ", r)
cex.cor <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = 2)
} ; upper.panel<-function(x, y){
  points(x,y, pch = 19)
}# Customize upper panel
pairs(GE3_3P_logFits.W[,2:8], #just another way to look at the data.
      lower.panel = upper.panel,
      upper.panel = panel.cor)
# GE3 Logistic GWAS ####
GE3Lower.mlmm = GAPIT(Y = as.data.frame(GE3_3P_logFits.W[,c('taxa','Lower')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                      Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3Lower.mlmm.s = DF_OperationsV2(GE3Lower.mlmm$GWAS)

GE3Rate.mlmm = GAPIT(Y = as.data.frame(GE3_3P_logFits.W[,c('taxa','Rate')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                     Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3Rate.mlmm.s = DF_OperationsV2(GE3Rate.mlmm$GWAS)

GE3Center.mlmm = GAPIT(Y = as.data.frame(GE3_3P_logFits.W[,c('taxa','Centering')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                       Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3Center.mlmm.s = DF_OperationsV2(GE3Center.mlmm$GWAS)

GE3Timeto90.mlmm = GAPIT(Y = as.data.frame(GE3_3P_logFits.W[,c('taxa','TimeTo90')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                         Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3Timeto90.mlmm.s = DF_OperationsV2(GE3Timeto90.mlmm$GWAS)

GE3Timeto95.mlmm = GAPIT(Y = as.data.frame(GE3_3P_logFits.W[,c('taxa','TimeTo95')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                         Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3Timeto95.mlmm.s = DF_OperationsV2(GE3Timeto95.mlmm$GWAS)

GE3rTimeto90.mlmm = GE3_3P_logFits.W %>% filter(rTimeTo90 < 200) %>% select(taxa,rTimeTo90)%>%data.frame()%>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3rTimeto90.mlmm.s = DF_OperationsV2(GE3rTimeto90.mlmm$GWAS)

GE3rTimeto95.mlmm = GE3_3P_logFits.W %>% filter(rTimeTo95 < 250) %>% select(taxa,rTimeTo95)%>%data.frame()%>%
  GAPIT(Y =.,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                          Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE3rTimeto95.mlmm.s = DF_OperationsV2(GE3rTimeto95.mlmm$GWAS)

GE3LogisticFits1920.GWAS.S = list(GE3Rate = GE3Rate.mlmm.s,
                                  GE3Lower = GE3Lower.mlmm.s,
                                  GE3Center=GE3Center.mlmm.s,
                                  GE3Timeto90 =GE3Timeto90.mlmm.s,
                                  GE3Timeto95= GE3Timeto95.mlmm.s,
                                  GE3rTimeto90= GE3rTimeto90.mlmm.s,
                                  GE3rTimeto95 =GE3rTimeto95.mlmm.s)
# get markers GE3 logistic GWAS #####
top_n_markers = function(GWAS.sum.df, trait,  ModelingType, n = 40){
  GWAS.sum.df %>% arrange(P.value) %>% filter(maf > 0.06) %>% slice_head(n = n) %>% 
    mutate(trait = trait,
           maf = round(as.numeric(as.character(maf)),2),
           ModelType = ModelingType)
}
GE3LogisticFits1920.GWAS.S.Traits = c('Rate','Lower','Centering','TimeTo90','TimeTo95','rTimeTo90','rTimeTo95')
GE3LogisticTopMarkers = data.frame()
counter = 1 
for (i in GE3LogisticFits1920.GWAS.S){
  GE3LogisticTopMarkers = rbind(GE3LogisticTopMarkers, 
                                top_n_markers(i,
                                              trait = GE3LogisticFits1920.GWAS.S.Traits[counter],
                                              ModelingType = 'GE3Logistic'))
  counter = counter + 1
}

SNPperChr = table(GE3LogisticFits1920.GWAS.S$GE3Rate$Chromosome)
MiddleOfChr = SNPperChr/2
breaks = c(); breaks[1]= MiddleOfChr[1];
ChromLines = c();ChromLines[1] = SNPperChr[1]
for(i in 2:length(SNPperChr)){
  breaks[i] = sum(SNPperChr[1:i-1])+MiddleOfChr[i]
  ChromLines[i] = sum(SNPperChr[1:i])
}

#  9228 is the ordinal position of the qds1 snp. 
#  ~10771 is the ordinal for the position of the Sd2 regions

GE3LogisticTopMarkers %>% 
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'GE3 Logistic MTA, MAF>0.06',
                  color = 'Parameter',shape = 'Model Type')

GE3LogisticTopMarkers %>% filter(P.value <5e-5 & maf >0.07) %>%
  arrange(trait,Chromosome, Position, P.value)  

GE3Marker_List = GE3LogisticTopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1)  

GE3Marker_List2 = GE3LogisticTopMarkers %>% filter(P.value <1e-4 & maf >0.07) %>%
  arrange(trait, Chromosome, Position) %>% select(SNP, Chromosome, Position, P.value, maf,trait)

myGD20_prune %>% select(GE3Marker_List$SNP) %>%
  ld_heatmap(.,GE3Marker_List)

# Repeat analysis with only the <8.0 line lower asymptopes #########
FiltGE3Lower.mlmm =GE3_3P_logFits.W  %>% filter(Lower<0.8) %>% select(taxa, Lower) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3Lower.mlmm.s = DF_OperationsV2(FiltGE3Lower.mlmm$GWAS)

FiltGE3Rate.mlmm = GE3_3P_logFits.W  %>% filter(Lower<0.8) %>% select(taxa, Rate) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3Rate.mlmm.s = DF_OperationsV2(FiltGE3Rate.mlmm$GWAS)

FiltGE3Center.mlmm = GE3_3P_logFits.W  %>% filter(Lower<0.8) %>% select(taxa, Centering) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3Center.mlmm.s = DF_OperationsV2(FiltGE3Center.mlmm$GWAS)

FiltGE3Timeto90.mlmm = GE3_3P_logFits.W  %>% filter(Lower<0.8) %>% select(taxa, TimeTo90) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3Timeto90.mlmm.s = DF_OperationsV2(GE3Timeto90.mlmm$GWAS)

FiltGE3Timeto95.mlmm = GE3_3P_logFits.W  %>% filter(Lower<0.8) %>% select(taxa, TimeTo95) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3Timeto95.mlmm.s = DF_OperationsV2(FiltGE3Timeto90.mlmm$GWAS)

FiltGE3rTimeto90.mlmm = GE3_3P_logFits.W %>% filter(rTimeTo90 < 200 & Lower < 0.8) %>% select(taxa,rTimeTo90)%>%data.frame()%>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3rTimeto90.mlmm.s = DF_OperationsV2(FiltGE3rTimeto90.mlmm$GWAS)

FiltGE3rTimeto95.mlmm = GE3_3P_logFits.W %>% filter(rTimeTo95 < 250 & Lower < 0.8) %>% select(taxa,rTimeTo95)%>%data.frame()%>%
  GAPIT(Y =.,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGE3rTimeto95.mlmm.s = DF_OperationsV2(FiltGE3rTimeto95.mlmm$GWAS)

FiltGE3LogisticFits1920.GWAS.S = list(fGE3Rate = FiltGE3Rate.mlmm.s,
                                  fGE3Lower = FiltGE3Lower.mlmm.s,
                                  fGE3Center=FiltGE3Center.mlmm.s,
                                  fGE3Timeto90 =FiltGE3Timeto90.mlmm.s,
                                  fGE3Timeto95= FiltGE3Timeto95.mlmm.s,
                                  fGE3rTimeto90= FiltGE3rTimeto90.mlmm.s,
                                  fGE3rTimeto95 =FiltGE3rTimeto95.mlmm.s)
filtGE3LogisticFits1920.GWAS.S.Traits = c('Rate','Lower','Centering','TimeTo90','TimeTo95','rTimeTo90','rTimeTo95')
filtGE3LogisticTopMarkers = data.frame()
counter = 1 
for (i in FiltGE3LogisticFits1920.GWAS.S){
  filtGE3LogisticTopMarkers = rbind(filtGE3LogisticTopMarkers, 
                                top_n_markers(i,
                                              trait = filtGE3LogisticFits1920.GWAS.S.Traits[counter],
                                              ModelingType = 'FiltGE3Logistic'))
  counter = counter + 1
}

filtGE3LogisticTopMarkers %>%  filter(maf >0.07 &P.value<5e-5) %>% 
  arrange(trait)%>%
  select(!c(nobs,effect, logPercentileQQplot, rank, FDRPval, ordinal))

# GE3 FPCA ##############
GE3.7Obs = All.BluesGE31920 %>% group_by(taxa) %>% dplyr::tally() %>% filter(n==7) %>% select(taxa)
All.BluesGE31920.FPCAInput = All.BluesGE31920 %>% filter(taxa %in% GE3.7Obs$taxa)
GE31920.FPCA = FPCA_function(dfTaxaTraitTime = All.BluesGE31920.FPCAInput[,c('taxa','time','GE3')],
                             Trait = 'GE3', #Trait name must be entered as a character ie Trait = 'GE3'
                             NumKnots = 1, # NumKnots is the number of interior knots to be fitted
                             order = 3, # Order is the dergree of the polynomial to be fit to the data.
                             NumObsevationsPerLine = 7)

GE31920.FPCA$EstimatedMeanPlot
GE31920.FPCA$PCs_withTaxa %>% select(FPC1, FPC2, FPC3) %>% pivot_longer(cols = starts_with('FPC'), names_to = 'FPC') %>%
  ggplot(aes(x = value))+geom_histogram()+facet_wrap(facets =vars(FPC),nrow = 1, scales = 'free')
GE31920.FPCA$RecoveredCurvePlot
  GE31920.FPCA$phi.fun.plot
  GE31920.FPCA$v1

VectorofTimeto95 = as.data.frame(GE31920.FPCA$RecoveredCurves) %>% mutate(time = GE31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:dim(GE31920.FPCA$PCs_withTaxa)[1],names_to = 'taxa') %>%
  group_by(taxa) %>% mutate(test  = value>0.95) %>% filter(test ==TRUE) %>% arrange(time) %>%slice_head()%>%ungroup() %>%
  mutate(TimeTo95 = time)

VectorofTimeto90 = as.data.frame(GE31920.FPCA$RecoveredCurves) %>% mutate(time = GE31920.FPCA$phi.fun.df$time) %>% 
  pivot_longer(cols = 1:dim(GE31920.FPCA$PCs_withTaxa)[1],names_to = 'taxa') %>%
  group_by(taxa) %>% mutate(test  = value>0.90) %>% filter(test ==TRUE) %>% arrange(time) %>%slice_head()%>%ungroup() %>%
  mutate(TimeTo90 = time)

GE31920FPCA.ouputs.longer = GE31920.FPCA$PCs_withTaxa %>% merge(.,VectorofTimeto90, by = 'taxa') %>% select(!c(time,value,test)) %>%
  merge(.,VectorofTimeto95, by = 'taxa') %>% select(!c(time,value,test)) %>%pivot_longer(cols = c(2:7)) 

GE31920FPCA.ouputs.longer %>% ggplot(aes(x = value))+geom_histogram()+facet_wrap(facets =vars(name), scales = 'free')

GE31920FPCA.ouputs.longer %>%pivot_wider(names_from = 'name',values_from = 'value') %>% select(!c(taxa,FPC4)) %>%
  pairs(., lower.panel = upper.panel, upper.panel = panel.cor)

GE31920PC1.mlmm = GAPIT(Y = GE31920.FPCA$PCs_withTaxa[,c('taxa','FPC1')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                             Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920PC1.mlmm.s = DF_OperationsV2(GE31920PC1.mlmm$GWAS)

GE31920PC2.mlmm = GAPIT(Y = GE31920.FPCA$PCs_withTaxa[,c('taxa','FPC2')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                             Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920PC2.mlmm.s = DF_OperationsV2(GE31920PC2.mlmm$GWAS)

GE31920PC3.mlmm = GAPIT(Y = GE31920.FPCA$PCs_withTaxa[,c('taxa','FPC3')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                             Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920PC3.mlmm.s = DF_OperationsV2(GE31920PC3.mlmm$GWAS)

GE31920FPCAtimeto95.mlmm = GAPIT(Y = data.frame(VectorofTimeto95[,c('taxa','TimeTo95')]) ,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                                      Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920FPCAtimeto95.mlmm.s = DF_OperationsV2(GE31920FPCAtimeto95.mlmm$GWAS)
TableOutput(GE31920FPCAtimeto95.mlmm.s)

GE31920FPCAtimeto90.mlmm = GAPIT(Y = data.frame(VectorofTimeto90[,c('taxa','TimeTo90')]) ,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                                      Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920FPCAtimeto90.mlmm.s = DF_OperationsV2(GE31920FPCAtimeto90.mlmm$GWAS)
TableOutput(GE31920FPCAtimeto90.mlmm.s)

GE3FPCA1920.GWAS.S = list(GE3FPCA_PC1 = GE31920PC1.mlmm.s,
                                  GE3FPCA_PC2 = GE31920PC2.mlmm.s,
                                  GE3FPCA_PC3 = GE31920PC3.mlmm.s,
                                  GE3FPCA_Timeto90 = GE31920FPCAtimeto90.mlmm.s,
                                  GE3FPCA_Timeto95 = GE31920FPCAtimeto95.mlmm.s)

GE3FPCA1920.GWAS.S.traits = c('FPC1','FPC2','FPC3','TimeTo90','TimeTo95')

GE3FPCATopMarkers = data.frame()
counter = 1 
for (i in GE3FPCA1920.GWAS.S){
  GE3FPCATopMarkers = rbind(GE3FPCATopMarkers, 
                                top_n_markers(i,
                                              trait = GE3FPCA1920.GWAS.S.traits[counter],
                                              ModelingType = 'GE3FPCA'))
  counter = counter + 1
}

GE3FPCATopMarkers %>% 
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'GE3 FPCA MTA, MAF>0.06',
                  color = 'Parameter',shape = 'Model Type')

GE3FPCAMarker_List = GE3FPCATopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1)  

myGD20_prune %>% select(GE3FPCAMarker_List$SNP) %>%
  ld_heatmap(.,GE3FPCAMarker_List)

# GI3 Logistic ###############
load('PhenotypeData/ProcessedData/2020/GGS2020_BLUE_summary_allTP.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM111_GI3_1920_BLUE.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM47_GI3_1920_BLUE.RData')
load('PhenotypeData/ProcessedData/2019and2020Combined/PM6_GI3_1920_BLUE.RData')

# Merge 2020 and 2019 data 
all1920GIBlues = PM6_GI3_1920_BLUE %>% mutate(TP = 'TP1') %>%
  rbind(.,(PM47_GI3_1920_BLUE %>% mutate(TP = 'TP4'))) %>%
  rbind(.,(PM111_GI3_1920_BLUE %>% mutate(TP = 'TP6'))) %>% 
  merge(., all_BLUE, by = c('taxa','TP'), all = TRUE) %>% 
  mutate(GI3 = ifelse(is.na(GI3),GI3scale, GI3),
         time = PM_date - 5) %>% 
  filter(!is.na(GE3)) %>%
  select(!c(GE3,GE5,GI3scale,GI5scale))#there are some errant lines that only appeared in 2019 data that must be removed. 
dim(all1920GIBlues)
all1920GIBlues %>% ggplot(aes(x = time, y = GI3, group = taxa))+geom_point()

# These make the loop break if all 7 timepoints are used
all1920GIBlues %>% filter(taxa %in%c('G_31_5', 'Megs_song','NY18120B_4','NY18125B_1','Oderbrucker',
                                     'P_2_5','P_26_3','P_36_6', 'SB193R_1','SG514R_4','SN873_3',
                                     'ST1431R_1')) 
all1920GIBlues %>% arrange(taxa, TP) %>% group_by(taxa) %>%  dplyr::tally() %>% select(n) %>% table() 

GI3GT5Obs = all1920GIBlues  %>% group_by(taxa) %>% dplyr::tally() %>% filter(n>4) %>% select(taxa) #GT = Greater than

GI3_4P_logFits = all1920GIBlues %>%
  filter(taxa %in% GI3GT5Obs$taxa & 
           taxa %nin% c('G_31_5', 'Megs_song','NY18120B_4','NY18125B_1','Oderbrucker',
                        'P_2_5','P_26_3','P_36_6', 'SB193R_1','SG514R_4','SN873_3',
                        'ST1431R_1')) %>%
  arrange(taxa, TP) %>% group_by(taxa) %>% 
  group_modify(~ broom::tidy( drm(GI3~time, data = .x,
                                  fct =LL.4(names = c('Rate','Lower','Upper','Centering'))))) %>%
  do(add_row(., taxa = .$taxa[1],term = 'TimeTo5.0',curve ='Derived',
             estimate = ifelse(.[2,4] > 5.0, 0,as.double((((.[3,4]-.[2,4])/(5.0-.[2,4])-1)^.[1,4])*.[4,4])))%>%
       add_row(., taxa = .$taxa[1],term = 'TimeTo5.6',curve ='Derived',
               estimate = ifelse(.[2,4] > 5.6, 0,as.double((((.[3,4]-.[2,4])/(5.6-.[2,4])-1)^.[1,4])*.[4,4])))%>%
       add_row(., taxa = .$taxa[1],term = 'DeltaGI95',curve ='Derived',
               estimate = as.double(.[4,4]*exp(log((1-0.95)/0.95)/.[1,4]))) %>%
       add_row(., taxa = .$taxa[1],term = 'DeltaGI90',curve ='Derived',
               estimate = as.double(.[4,4]*exp(log((1-0.90)/0.90)/.[1,4])))
     ) %>% 
  ungroup() %>% mutate(estimate = ifelse(is.nan(estimate),NA,estimate))
# We give any models that fail to reach a GI of that level a NA at this point. 
GI3_4P_logFits %>% ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')
GI3_4P_logFits %>% filter(!(term == 'Lower' & estimate > 10 |
                            term == 'Rate' &(estimate > 0 | estimate < -15) |
                            term == 'Centering' & estimate > 150)) %>% 
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')

GI3logTaxatoFilter = GI3_4P_logFits %>% 
  filter((term == 'Lower' & estimate > 10 |
            term == 'Rate' &(estimate > 0 | estimate < -15) |
            term == 'Centering' & estimate > 150)) %>% select(taxa) %>% unique()

GI3fittedCurves = data.frame()
for (i in unique(GI3_4P_logFits$taxa)){
  time = seq(0,150,3)
  tmp = GI3_4P_logFits %>% filter(taxa ==i)
  y = tmp$estimate[2]+(tmp$estimate[3]-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[4]))))
  GI3fittedCurves = rbind(GI3fittedCurves, data.frame(taxa = i, time = time, GI3estimate = y))
}
GI3fittedCurves %>% ggplot(aes(x = time, y = GI3estimate, group = taxa))+
  geom_line() +labs(title = 'GI3 Logistic Fits')

GI3fittedCurves %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% ggplot(aes(x = time, y = GI3estimate, group = taxa))+
  geom_line() +labs(title = 'GI3 Logistic Fits')
GI3_4P_logFits %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% 
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')
GI3_4P_logFits %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% filter(term %in% c('Centering','Upper','Lower','Rate')) %>%
  ggplot(aes( x= estimate))+geom_histogram()+facet_wrap(vars(term), scales = 'free')

GI3_4P_logFits.w = GI3_4P_logFits %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% 
  select(taxa, term, estimate) %>% pivot_wider(names_from = 'term', values_from = 'estimate')%>%
  mutate(DeltaGI = Upper-Lower,
         DeltaGI90 = ifelse(DeltaGI90>150,NA,DeltaGI90),
         DeltaGI95 = ifelse(DeltaGI95>200,NA,DeltaGI95)) 

GI3Lower.mlmm = GAPIT(Y = as.data.frame(GI3_4P_logFits.w[,c('taxa','Lower')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                      Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3Lower.mlmm.s = DF_OperationsV2(GI3Lower.mlmm$GWAS)

GI3Rate.mlmm = GAPIT(Y = as.data.frame(GI3_4P_logFits.w[,c('taxa','Rate')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                     Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3Rate.mlmm.s = DF_OperationsV2(GI3Rate.mlmm$GWAS)

GI3Center.mlmm = GAPIT(Y = as.data.frame(GI3_4P_logFits.w[,c('taxa','Centering')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                       Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3Center.mlmm.s = DF_OperationsV2(GI3Center.mlmm$GWAS)

GI3Upper.mlmm = GAPIT(Y = as.data.frame(GI3_4P_logFits.w[,c('taxa','Upper')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                       Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3Upper.mlmm.s = DF_OperationsV2(GI3Upper.mlmm$GWAS)

GI3DeltaGI.mlmm = GAPIT(Y = as.data.frame(GI3_4P_logFits.w[,c('taxa','DeltaGI')]),GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,
                      Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3DeltaGI.mlmm.s = DF_OperationsV2(GI3DeltaGI.mlmm$GWAS)

GI3DeltaGI90.mlmm = GI3_4P_logFits.w %>% select(taxa, DeltaGI90) %>% filter(DeltaGI90<150)%>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3DeltaGI90.mlmm.s = DF_OperationsV2(GI3DeltaGI90.mlmm$GWAS)
TableOutput(GI3DeltaGI90.mlmm.s)

GI3DeltaGI95.mlmm = GI3_4P_logFits.w %>% select(taxa, DeltaGI95) %>% filter(DeltaGI95<200)%>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3DeltaGI95.mlmm.s = DF_OperationsV2(GI3DeltaGI95.mlmm$GWAS)
TableOutput(GI3DeltaGI95.mlmm.s)

GI3Timeto5.0.mlmm = GI3_4P_logFits.w %>% select(taxa, TimeTo5.0) %>% filter(TimeTo5.0<250) %>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI3Timeto5.0.mlmm.s = DF_OperationsV2(GI3Timeto5.0.mlmm$GWAS)
TableOutput(GI3Timeto5.0.mlmm.s)

GI3LogisticFits1920.GWAS.S = list(GI3Rate = GI3Rate.mlmm.s,
                                  GI3Lower = GI3Lower.mlmm.s,
                                  GI3Center=GI3Center.mlmm.s,
                                  GI3Upper = GI3Upper.mlmm.s,
                                  GI3DeltaGI= GI3DeltaGI.mlmm.s,
                                  GI3DeltaGI90 =GI3DeltaGI90.mlmm.s,
                                  GI3DeltaGI95 = GI3DeltaGI95.mlmm.s,
                                  GI3Timeto5.0 = GI3Timeto5.0.mlmm.s)
GI3LogisticFits1920.GWAS.S.Traits = c('Rate','Lower','Centering','Upper',
                                      'DeltaGI','DeltaGI90','DeltaGI95', 'TimeTo5.0')

GI3LogisticTopMarkers = data.frame()
counter = 1 
for (i in GI3LogisticFits1920.GWAS.S){
  GI3LogisticTopMarkers = rbind(GI3LogisticTopMarkers, 
                                top_n_markers(i,
                                              trait = GI3LogisticFits1920.GWAS.S.Traits[counter],
                                              ModelingType = 'GI3Logistic'))
  counter = counter + 1
}

GI3LogisticTopMarkers %>% 
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'GI3 Logistic MTA, MAF>0.06',
                  color = 'Parameter',shape = 'Model Type')
GI3LogisticTopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1 | row_number()==n())  

GI3Marker_List = GI3LogisticTopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1)  

myGD20_prune %>% select(GI3Marker_List$SNP) %>%
  ld_heatmap(.,GI3Marker_List)
 
pairs(GI3_4P_logFits.w[,2:5], 
      lower.panel = upper.panel,
      upper.panel = panel.cor)
# GI3 GWA on parameters after filtering for GI>5 #####
FiltGI3Lower.mlmm =GI3_4P_logFits.w  %>% filter(Lower<5) %>% select(taxa, Lower) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3Lower.mlmm.s = DF_OperationsV2(FiltGI3Lower.mlmm$GWAS)

FiltGI3Rate.mlmm = GI3_4P_logFits.w  %>% filter(Lower<5) %>% select(taxa, Rate) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3Rate.mlmm.s = DF_OperationsV2(FiltGI3Rate.mlmm$GWAS)

FiltGI3Center.mlmm = GI3_4P_logFits.w  %>% filter(Lower<5) %>% select(taxa, Centering) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3Center.mlmm.s = DF_OperationsV2(FiltGI3Center.mlmm$GWAS)

FiltGI3Upper.mlmm = GI3_4P_logFits.w  %>% filter(Lower<5) %>% select(taxa, Upper) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3Upper.mlmm.s = DF_OperationsV2(FiltGI3Upper.mlmm$GWAS)

FiltGI3DeltaGI.mlmm = GI3_4P_logFits.w  %>% filter(Lower<5) %>% select(taxa, DeltaGI) %>% as.data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3DeltaGI.mlmm.s = DF_OperationsV2(FiltGI3DeltaGI.mlmm$GWAS)

FiltGI3DeltaGI90.mlmm = GI3_4P_logFits.w %>% filter(Lower<5) %>% select(taxa, DeltaGI90) %>% filter(DeltaGI90<150)%>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3DeltaGI90.mlmm.s = DF_OperationsV2(FiltGI3DeltaGI90.mlmm$GWAS)

FiltGI3DeltaGI95.mlmm = GI3_4P_logFits.w %>% filter(Lower<5) %>% select(taxa, DeltaGI95) %>% filter(DeltaGI95<200)%>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3DeltaGI95.mlmm.s = DF_OperationsV2(FiltGI3DeltaGI95.mlmm$GWAS)

FiltGI3Timeto5.0.mlmm = GI3_4P_logFits.w %>% filter(Lower<5) %>% select(taxa, TimeTo5.0) %>% filter(TimeTo5.0<250) %>% data.frame() %>%
  GAPIT(Y = .,GD=myGD20_prune, GM=myGM20_prune,PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
FiltGI3Timeto5.0.mlmm.s = DF_OperationsV2(FiltGI3Timeto5.0.mlmm$GWAS)

FiltGI3LogisticFits1920.GWAS.S = list(fGI3Rate = FiltGI3Rate.mlmm.s,
                                      fGI3Lower = FiltGI3Lower.mlmm.s,
                                      fGI3Center=FiltGI3Center.mlmm.s,
                                      fGI3Upper = FiltGI3Upper.mlmm.s,
                                      fGI3DeltaGI = FiltGI3DeltaGI.mlmm.s,
                                      fGI3DetlaGI90 = FiltGI3DeltaGI90.mlmm.s,
                                      fGI3DetlaGI95 = FiltGI3DeltaGI95.mlmm.s,
                                      fGI3TimeTo5.0 = FiltGI3Timeto5.0.mlmm.s)
filtGI3LogisticFits1920.GWAS.S.Traits = c('Rate','Lower','Centering','Upper',
                                          'DeltaGI','DeltaGI90','DeltaGI95','TimeTo5.0')
filtGI3LogisticTopMarkers = data.frame()
counter = 1 
for (i in FiltGI3LogisticFits1920.GWAS.S){
  filtGI3LogisticTopMarkers = rbind(filtGI3LogisticTopMarkers, 
                                    top_n_markers(i,
                                                  trait = filtGI3LogisticFits1920.GWAS.S.Traits[counter],
                                                  ModelingType = 'FiltGI3Logistic'))
  counter = counter + 1
}

filtGI3LogisticTopMarkers %>%  filter(maf >0.07 &P.value<5e-5) %>% 
  arrange(trait)%>%
  select(!c(nobs,effect, logPercentileQQplot, rank, FDRPval, ordinal))


# GI3 FPCA #################################################

all1920GIBlues$time
GI3.7Obs = all1920GIBlues %>% group_by(taxa) %>% dplyr::tally() %>% filter(n==7) %>% select(taxa)
all1920GIBlues.FPCAInput = all1920GIBlues %>% filter(taxa %in% GI3.7Obs$taxa)
GI31920.FPCA = FPCA_function(dfTaxaTraitTime = all1920GIBlues.FPCAInput[,c('taxa','time','GI3')],
                             Trait = 'GI3', #Trait name must be entered as a character ie Trait = 'GE3'
                             NumKnots = 1, # NumKnots is the number of interior knots to be fitted
                             order = 3, # Order is the dergree of the polynomial to be fit to the data.
                             NumObsevationsPerLine = 7)

GI31920.FPCA$EstimatedMeanPlot
GI31920.FPCA$PCs_withTaxa %>% select(FPC1, FPC2, FPC3) %>% pivot_longer(cols = starts_with('FPC'), names_to = 'FPC') %>%
  ggplot(aes(x = value))+geom_histogram()+facet_wrap(facets =vars(FPC),nrow = 1, scales = 'free')
GI31920.FPCA$RecoveredCurvePlot

VectorofTimeto5.0 = as.data.frame(GI31920.FPCA$RecoveredCurves) %>% mutate(time = GI31920.FPCA$phi.fun.df$time) %>% 
  pivot_longer(cols = 1:dim(GI31920.FPCA$PCs_withTaxa)[1],names_to = 'taxa') %>%
  group_by(taxa) %>% mutate(test  = value>5.0) %>% filter(test ==TRUE) %>% arrange(time) %>%slice_head()%>%ungroup() %>%
  mutate(TimeTo5.0 = time)
hist(VectorofTimeto5.0$TimeTo5.0)

VectorofTimeto5.6 = as.data.frame(GI31920.FPCA$RecoveredCurves) %>% mutate(time = GI31920.FPCA$phi.fun.df$time) %>% 
  pivot_longer(cols = 1:dim(GI31920.FPCA$PCs_withTaxa)[1],names_to = 'taxa') %>%
  group_by(taxa) %>% mutate(test  = value>5.6) %>% filter(test ==TRUE) %>% arrange(time) %>%slice_head()%>%ungroup() %>%
  mutate(TimeTo5.6 = time)
hist(VectorofTimeto5.6$TimeTo5.6)
GI31920FPCA.ouputs.longer = GI31920.FPCA$PCs_withTaxa %>% merge(.,VectorofTimeto5.0, by = 'taxa',all.x = TRUE) %>% select(!c(time,value,test)) %>%
  merge(.,VectorofTimeto5.6, by = 'taxa',all.x = TRUE) %>% select(!c(time,value,test)) %>%pivot_longer(cols = c(2:7)) 

GI31920FPCA.ouputs.longer %>% ggplot(aes(x = value))+
  geom_histogram()+facet_wrap(facets =vars(name), scales = 'free')

GI31920FPCA.ouputs.longer %>%pivot_wider(names_from = 'name',values_from = 'value') %>% select(!c(taxa,FPC4)) %>%
  pairs(., 
        lower.panel = upper.panel,
        upper.panel = panel.cor)

GI31920.FPCA$v1

GI31920PC1.mlmm = GAPIT(Y = GI31920.FPCA$PCs_withTaxa[,c('taxa','FPC1')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                        Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920PC1.mlmm.s = DF_OperationsV2(GI31920PC1.mlmm$GWAS)

GI31920PC2.mlmm = GAPIT(Y = GI31920.FPCA$PCs_withTaxa[,c('taxa','FPC2')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                        Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920PC2.mlmm.s = DF_OperationsV2(GI31920PC2.mlmm$GWAS)

GI31920PC3.mlmm = GAPIT(Y = GI31920.FPCA$PCs_withTaxa[,c('taxa','FPC3')],GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                        Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920PC3.mlmm.s = DF_OperationsV2(GI31920PC3.mlmm$GWAS)

GI31920FPCAtimeto5.0.mlmm = GAPIT(Y = data.frame(VectorofTimeto5.0[,c('taxa','TimeTo5.0')]) ,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                                 Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920FPCAtimeto5.0.mlmm.s = DF_OperationsV2(GI31920FPCAtimeto5.0.mlmm$GWAS)
TableOutput(GI31920FPCAtimeto5.0.mlmm.s)

GI31920FPCAtimeto5.6.mlmm = GAPIT(Y = data.frame(VectorofTimeto5.6[,c('taxa','TimeTo5.6')]) ,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,
                                 Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920FPCAtimeto5.6.mlmm.s = DF_OperationsV2(GI31920FPCAtimeto5.6.mlmm$GWAS)
TableOutput(GI31920FPCAtimeto5.6.mlmm.s)

GI3FPCA1920.GWAS.S = list(GI3FPCA_PC1 = GI31920PC1.mlmm.s,
                          GI3FPCA_PC2 = GI31920PC2.mlmm.s,
                          GI3FPCA_PC3 = GI31920PC3.mlmm.s,
                          GI3FPCA_Timeto5.0 = GI31920FPCAtimeto5.0.mlmm.s,
                          GI3FPCA_Timeto5.6 = GI31920FPCAtimeto5.6.mlmm.s)

GI3FPCA1920.GWAS.S.traits = c('FPC1','FPC2','FPC3','TimeTo5.0','TimeTo5.6')

GI3FPCATopMarkers = data.frame()
counter = 1 
for (i in GI3FPCA1920.GWAS.S){
  GI3FPCATopMarkers = rbind(GI3FPCATopMarkers, 
                            top_n_markers(i,
                                          trait = GI3FPCA1920.GWAS.S.traits[counter],
                                          ModelingType = 'GI3FPCA'))
  counter = counter + 1
}

GI3FPCATopMarkers %>% 
  ggplot(aes(x = ordinal, y = log10PVal))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'GI3 FPCA MTA, MAF>0.06',
                  color = 'Parameter',shape = 'Model Type')

GI3FPCATopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1 | row_number()==n())  

GI3FPCAMarker_List = GI3FPCATopMarkers %>% filter(P.value <1e-4) %>%
  arrange(Chromosome, Position, P.value) %>% group_by(SNP) %>%
  filter(row_number()==1)  

myGD20_prune %>% select(GI3FPCAMarker_List$SNP) %>%
  ld_heatmap(.,GI3FPCAMarker_List)


# Correlation between TP1 to 7 values for GI and GE and the time series values predicted for those times #####
# 0  14  28  42  63 105 154 ARE THE TIMES THAT things were measured,so we need close to that...
GE31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = round(GE31920.FPCA$phi.fun.df$time,2)) %>%
  pivot_longer(cols = 1:484, values_to = 'GE3', names_to = 'taxa') %>% select(taxa,time,GE3) %>%
  filter(time %in% c(0,13.91,27.81,42.23,62.84,105.07,154)) %>% 
  mutate(time = mapvalues(time, from =c(0,13.91,27.81,42.23,62.84,105.07,154),to = c(0,14,28,42,63,105,154))) %>% rename(fGE = GE3) %>%
  merge(All.BluesGE31920.FPCAInput, by=c('taxa','time')) %>% group_by(time) %>%summarise(correlation = cor(fGE,GE3))

GI31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = round(GI31920.FPCA$phi.fun.df$time,2)) %>%
  pivot_longer(cols = 1:484, values_to = 'GI3', names_to = 'taxa') %>% select(taxa,time,GI3) %>%
  filter(time %in% c(0,13.91,27.81,42.23,62.84,105.07,154)) %>% 
  mutate(time = mapvalues(time, from =c(0,13.91,27.81,42.23,62.84,105.07,154),to = c(0,14,28,42,63,105,154))) %>% rename(fGI = GI3) %>%
  merge(all1920GIBlues.FPCAInput, by=c('taxa','time')) %>% group_by(time) %>%summarise(correlation = cor(fGI,GI3))

GE3logEstimatesForComparison = data.frame()
for (i in unique(GE3_3P_logFits$taxa)){
  time =c(0,14,28,42,63,105,154)
  tmp = GE3_3P_logFits %>% filter(taxa ==i)
  y = tmp$estimate[2]+(1-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[3]))))
  GE3logEstimatesForComparison = rbind(GE3logEstimatesForComparison, data.frame(taxa = i, time = time, GE3estimate = y))
}
GE3logEstimatesForComparison %>% merge(.,All.BluesGE31920 %>% 
                                         filter(taxa %in% GE3GT4Obs$taxa & 
                                                  taxa %nin% c('Chevallier_D10',"G_31_5","P_16_1",'P_29_5')), by = c('taxa','time')) %>%
  group_by(time) %>% summarise(correlation = cor(GE3, GE3estimate))

GI3logEstimatesForComparison = data.frame()
for (i in unique(GE3_3P_logFits$taxa)){
  time =c(0,14,28,42,63,105,154)
  tmp = GI3_4P_logFits %>% filter(taxa ==i)
  y = tmp$estimate[2]+(tmp$estimate[3]-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[4]))))
  GI3logEstimatesForComparison = rbind(GI3logEstimatesForComparison, data.frame(taxa = i, time = time, GI3estimate = y))
}
GI3logEstimatesForComparison %>%  merge(.,all1920GIBlues %>% 
                                          filter(taxa %in% GI3GT5Obs$taxa & 
                                                   taxa %nin% c('G_31_5', 'Megs_song','NY18120B_4','NY18125B_1','Oderbrucker',
                                                                'P_2_5','P_26_3','P_36_6', 'SB193R_1','SG514R_4','SN873_3','ST1431R_1')),
                                        by = c('taxa','time')) %>%
  group_by(time) %>% summarise(correlation = cor(GI3, GI3estimate))

  
  
  
  
  


# Sum total of markers  #########
rbind(GI3FPCATopMarkers, GE3FPCATopMarkers,
      GI3LogisticTopMarkers,
      GE3LogisticTopMarkers) %>%  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait)) %>%
  ggplot(aes(x = ordinal, y = log10PVal, alpha = maf))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'All MTA, MAF>0.07',
                  color = 'Parameter',shape = 'Model Type')


SigMarkers = rbind(GI3FPCATopMarkers, GE3FPCATopMarkers,
                   GI3LogisticTopMarkers,
                   GE3LogisticTopMarkers) %>%  filter(maf >0.1 & P.value < 5e-5) %>%
  arrange(Chromosome, Position)
SigMarkers %>% arrange(Chromosome, Position, P.value) %>%select(!c(nobs,effect, logPercentileQQplot, rank, FDRPval, ordinal))
u = SigMarkers %>% select(SNP) %>% unique()

# Single Time point GWAS for Comparison #####
# All.BluesGE31920 is input
GE31920.TP1.mlmm = All.BluesGE31920 %>% filter(TP =='TP1') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP1.mlmm.s = DF_OperationsV2(GE31920.TP1.mlmm$GWAS)

GE31920.TP2.mlmm = All.BluesGE31920 %>% filter(TP =='TP2') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP2.mlmm.s = DF_OperationsV2(GE31920.TP2.mlmm$GWAS)

GE31920.TP3.mlmm = All.BluesGE31920 %>% filter(TP =='TP3') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP3.mlmm.s = DF_OperationsV2(GE31920.TP3.mlmm$GWAS)

GE31920.TP4.mlmm = All.BluesGE31920 %>% filter(TP =='TP4') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP4.mlmm.s = DF_OperationsV2(GE31920.TP4.mlmm$GWAS)

GE31920.TP5.mlmm = All.BluesGE31920 %>% filter(TP =='TP5') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP5.mlmm.s = DF_OperationsV2(GE31920.TP5.mlmm$GWAS)

GE31920.TP6.mlmm = All.BluesGE31920 %>% filter(TP =='TP6') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP6.mlmm.s = DF_OperationsV2(GE31920.TP6.mlmm$GWAS)

GE31920.TP7.mlmm = All.BluesGE31920 %>% filter(TP =='TP7') %>% select(taxa, GE3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GE31920.TP7.mlmm.s = DF_OperationsV2(GE31920.TP7.mlmm$GWAS)

GE31920PerTP.GWAS.S = list(
  TP1GE31920 = GE31920.TP1.mlmm.s,
  TP2GE31920 = GE31920.TP2.mlmm.s,
  TP3GE31920 = GE31920.TP3.mlmm.s,
  TP4GE31920 = GE31920.TP4.mlmm.s,
  TP5GE31920 = GE31920.TP5.mlmm.s,
  TP6GE31920 = GE31920.TP6.mlmm.s,
  TP7GE31920 = GE31920.TP7.mlmm.s
)
GE3perTP.GWAS.traits = c('TP1 GE3','TP2 GE3','TP3 GE3','TP4 GE3','TP5 GE3','TP6 GE3','TP7 GE3')
GE3perTPTopMarkers = data.frame()
counter = 1 
for (i in GE31920PerTP.GWAS.S){
  GE3perTPTopMarkers = rbind(GE3perTPTopMarkers, 
                            top_n_markers(i,
                                          trait = GE3perTP.GWAS.traits[counter],
                                          ModelingType = 'GE3 per Time Point'))
  counter = counter + 1
}

GE3perTPTopMarkers %>%  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait)) %>%
  ggplot(aes(x = ordinal, y = log10PVal, alpha = maf))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'All MTA, MAF>0.07',
                  color = 'Parameter',shape = 'Model Type')

# all1920GIBlues GI3 is input
GI31920.TP1.mlmm = all1920GIBlues %>% filter(TP =='TP1') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP1.mlmm.s = DF_OperationsV2(GI31920.TP1.mlmm$GWAS)

GI31920.TP2.mlmm = all1920GIBlues %>% filter(TP =='TP2') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP2.mlmm.s = DF_OperationsV2(GI31920.TP2.mlmm$GWAS)

GI31920.TP3.mlmm = all1920GIBlues %>% filter(TP =='TP3') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP3.mlmm.s = DF_OperationsV2(GI31920.TP3.mlmm$GWAS)

GI31920.TP4.mlmm = all1920GIBlues %>% filter(TP =='TP4') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP4.mlmm.s = DF_OperationsV2(GI31920.TP4.mlmm$GWAS)

GI31920.TP5.mlmm = all1920GIBlues %>% filter(TP =='TP5') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP5.mlmm.s = DF_OperationsV2(GI31920.TP5.mlmm$GWAS)

GI31920.TP6.mlmm = all1920GIBlues %>% filter(TP =='TP6') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP6.mlmm.s = DF_OperationsV2(GI31920.TP6.mlmm$GWAS)

GI31920.TP7.mlmm = all1920GIBlues %>% filter(TP =='TP7') %>% select(taxa, GI3) %>%data.frame() %>%
  GAPIT(.,GD=myGD20_prune, GM=myGM20_prune, PCA.total = 2,Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
GI31920.TP7.mlmm.s = DF_OperationsV2(GI31920.TP7.mlmm$GWAS)

GI31920PerTP.GWAS.S = list(
  TP1GI31920 = GI31920.TP1.mlmm.s,
  TP2GI31920 = GI31920.TP2.mlmm.s,
  TP3GI31920 = GI31920.TP3.mlmm.s,
  TP4GI31920 = GI31920.TP4.mlmm.s,
  TP5GI31920 = GI31920.TP5.mlmm.s,
  TP6GI31920 = GI31920.TP6.mlmm.s,
  TP7GI31920 = GI31920.TP7.mlmm.s
)
GI3perTP.GWAS.traits = c('TP1 GI3','TP2 GI3','TP3 GI3','TP4 GI3','TP5 GI3','TP6 GI3','TP7 GI3')
GI3perTPTopMarkers = data.frame()
counter = 1 
for (i in GI31920PerTP.GWAS.S){
  GI3perTPTopMarkers = rbind(GI3perTPTopMarkers, 
                             top_n_markers(i,
                                           trait = GI3perTP.GWAS.traits[counter],
                                           ModelingType = 'GI3 per Time Point'))
  counter = counter + 1
}

GI3perTPTopMarkers %>%  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait)) %>%
  ggplot(aes(x = ordinal, y = log10PVal, alpha = maf))+
  geom_point(aes(shape = ModelType, color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p)')+xlab(NULL)+
  geom_hline(yintercept = 4)+
  ylim(0,30)+
  theme_bw()+labs(title = 'All MTA, MAF>0.07',
                  color = 'Parameter',shape = 'Model Type')


setwd(rprojroot::find_rstudio_root_file())
# To save all the outputs uncomment and change the dirrectory to where you would like it
# setwd('TimeSeriesAnalysis/Output/')
# save(GI3FPCA1920.GWAS.S, file = 'GI3FPCA1920.GWAS.S.RData')
# save(GE3FPCA1920.GWAS.S, file = 'GE3FPCA1920.GWAS.S.RData')
# save(GI3LogisticFits1920.GWAS.S, file = 'GI3LogisticFits1920.GWAS.S.RData')
# save(GE3LogisticFits1920.GWAS.S, file = 'GE3LogisticFits1920.GWAS.S.RData')
# save(GI31920PerTP.GWAS.S, file = 'GI31920PerTP.GWAS.S.RData')
# save(GE31920PerTP.GWAS.S, file = 'GE31920PerTP.GWAS.S.RData')
# save(FiltGE3LogisticFits1920.GWAS.S, file ='FiltGE3LogisticFits1920.GWAS.S.RData')
# save(FiltGI3LogisticFits1920.GWAS.S, file = 'FiltGI3LogisticFits1920.GWAS.S.RData')
# setwd(rprojroot::find_rstudio_root_file())

##############################
### Genomic Prediction #######
setwd(rprojroot::find_rstudio_root_file())
load('SpringBarley/GenotypeData/myGD_LDpruned_w_KASP.RData')
load('SpringBarley/PhenotypeData/ProcessedData/2020/GGS2020_BLUE_summary_allTP.RData')
# We are interested in the GP power within our 
# breeding program so the JIC and naked barleys were excluded
PopulationKey = read.csv(file = 'SpringBarley/GenotypeData/PopulationKey.csv') %>% select(!X)
CULines =PopulationKey %>% filter(Population %in% c("check","C1G","C2G","C1P","base"))
CUTP1Germs = all_BLUE %>% filter(TP == "TP1" & taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6'))

myGD20_pruneCU = myGD20_prune %>% filter(taxa %in% CULines$taxa)
MAFList = data.frame()
counter = 2 
MarkerNames = colnames(myGD20_pruneCU)
for (col in myGD20_pruneCU[,-1]){
  MAFList = rbind(MAFList, data.frame(MarkerName = MarkerNames[counter], 
                                      maf = sum(as.numeric(col))/422/2))
  counter = counter+1
}
MAFList =MAFList %>% mutate(filt95 = maf<0.95, filt90 = maf<0.9)


myGD20_pruneCUfilt = myGD20_pruneCU[,which(c(TRUE,MAFList$filt95))]
myGD20_pruneCUfilt[1:5,1:5]
#takes 91 out to drop mkk3 
Mkk3taxa = myGD20_pruneCUfilt %>% select(taxa, MKK3_E165Q) %>% filter(MKK3_E165Q==0) %>% select(taxa)
dim(myGD20_pruneCUfilt)
myGM20_pruneCU = myGM20_prune[which(MAFList$filt95),]

FoldCVGP = function(Phenotype, myGD, numfolds=50,datasplit, trait_name){
  # Phenotype is full vector of phenotypes
  # myGD is a n taxa x N markers dataframe with -1,0,1 coding and no 'taxa' column in it
  # numfolds is the number of folds to cv (ie 5 for five-fold cv)
  # datasplit is the percentage as a decimal for the training/testing split.
  
  TrueandPredicted = data.frame()
  num_entries = length(Phenotype)
  Training_size = round(datasplit*num_entries,0)
  start_time <- Sys.time()
  set.seed(1)
  for (i in 1:numfolds){
    
    trainingSet = sort(sample(1:num_entries,Training_size))
    testingSet = setdiff(1:num_entries,trainingSet)
    
    y_train = Phenotype[trainingSet] 
    y_test = Phenotype[testingSet]
    
    marker_train = myGD[trainingSet,]
    marker_test = myGD[testingSet,]
    
    trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
    
    PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)
    print(cor(y_test,PredictedPheno, use = 'complete.obs'))
    
    TrueandPredicted = rbind(TrueandPredicted, 
                             data.frame(TruePheno = y_test,
                                        PredPheno = PredictedPheno,
                                        fold = i,
                                        trait = trait_name,
                                        correlation = cor(y_test,PredictedPheno, use = 'complete.obs')))
  }
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}

FoldCVGPSamePopSize = function(Phenotype1, myGD1, numfolds=50,datasplit, trait_name, sample_size){
  # Phenotype is full vector of phenotypes
  # myGD is a n taxa x N markers dataframe with -1,0,1 coding and no 'taxa' column in it
  # numfolds is the number of folds to cv (ie 5 for five-fold cv)
  # datasplit is the percentage as a decimal for the training/testing split.
  
  set.seed(1)
  TrueandPredicted = data.frame()
  num_entries = length(Phenotype1)
  Training_size = round(datasplit*sample_size,0)
  start_time <- Sys.time()
  
  for (i in 1:numfolds){
    subsampled = sort(sample(1:num_entries, sample_size))
    Phenotype = Phenotype1[subsampled]
    myGD = myGD1[subsampled,]
    
    trainingSet = sort(sample(1:sample_size,Training_size))
    testingSet = setdiff(1:sample_size,trainingSet)
    
    y_train = Phenotype[trainingSet] 
    y_test = Phenotype[testingSet]
    
    marker_train = myGD[trainingSet,]
    marker_test = myGD[testingSet,]
    
    trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
    
    PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)
    print(cor(y_test,PredictedPheno, use = 'complete.obs'))
    
    TrueandPredicted = rbind(TrueandPredicted, 
                             data.frame(TruePheno = y_test,
                                        PredPheno = PredictedPheno,
                                        fold = i,
                                        trait = trait_name,
                                        correlation = cor(y_test,PredictedPheno, use = 'complete.obs')))
  }
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}



GPPHSCU = rbind(FoldCVGP((PHS.blues %>% filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6')))$PHS,
                         myGD = myGD20_pruneCUfilt[,-1]-1, #convert to correct format
                         datasplit = .8,
                         trait_name = 'perTP_PHS_withVND'),
                FoldCVGP((PHS.blues %>% filter(taxa %in% CULines$taxa & 
                                                 taxa %nin% c('P_5_3','P_8_6') &
                                                 taxa %in% Mkk3taxa$taxa))$PHS,
                         myGD = (myGD20_pruneCUfilt %>% filter(taxa %in% Mkk3taxa$taxa))[,-1]-1, #convert to correct format
                         datasplit = .8,
                         trait_name = 'perTP_PHS_noVND'),
                FoldCVGPSamePopSize((PHS.blues %>% filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6')))$PHS,
                                    myGD = myGD20_pruneCUfilt[,-1]-1, #convert to correct format
                                    datasplit = .8,
                                    trait_name = 'perTP_PHS_withVNDSamePop',
                                    sample_size = length((PHS.blues %>% filter(taxa %in% CULines$taxa & 
                                                                                 taxa %nin% c('P_5_3','P_8_6') &
                                                                                 taxa %in% Mkk3taxa$taxa))$PHS)))

GPPHSCU %>% group_by(trait) %>% summarise(mean(correlation))
# GP Per TP GP using all the 1920BLUES.#####
# All.BluesGE31920 # contians the GE3
# all1920GIBlues # contains the GI3
all_BLUE.w =All.BluesGE31920 %>% select(!time) %>% mutate(trait = 'GE3') %>%
  pivot_wider(names_from =c( 'TP', 'trait'), values_from = 'GE3') %>%
  join(all1920GIBlues %>% select(!c('PM_date','time','date')) %>% mutate(trait = 'GI3scale') %>%
         pivot_wider(names_from =c( 'TP', 'trait'), values_from = 'GI3' ))  %>%
  filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6'))
PhenoNames = gsub(gsub(x = colnames(all_BLUE.w),pattern = '_', replacement = ''),pattern = 'scale',replacement = '')
GPPerTP = data.frame()

counter = 2
for (col in all_BLUE.w[,-1]){
  temp =  FoldCVGP(col,
                   myGD20_pruneCUfilt[,-1]-1, #convert to correct format
                   datasplit = .8,
                   trait_name = paste0('perTP_',PhenoNames[counter],'_withVND'))
  temp2 = FoldCVGP(col[which(all_BLUE.w$taxa %in% Mkk3taxa$taxa)],
                   (myGD20_pruneCUfilt %>% filter(taxa %in% Mkk3taxa$taxa &
                                                    taxa %in% all_BLUE.w$taxa))[,-1]-1,
                   datasplit = .8,
                   trait_name = paste0('perTP_',PhenoNames[counter],'_noVND'))
  temp3 = FoldCVGPSamePopSize(col,
                              myGD20_pruneCUfilt[,-1]-1, #convert to correct format
                              datasplit = .8,
                              trait_name = paste0('perTP_',PhenoNames[counter],'_withVNDSamePop'),
                              sample_size = length(col[which(all_BLUE.w$taxa %in% Mkk3taxa$taxa)]))
  GPPerTP = rbind(GPPerTP,temp,temp2,temp3)
  counter = counter+1
}

# GP GE3 Logistic fits #####
GPGElogfits = data.frame()
GE3_3P_logFits.WCU = GE3_3P_logFits.W %>% filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6'))
missingGE3Plog = myGD20_pruneCUfilt$taxa %in% GE3_3P_logFits.WCU$taxa
counter = 2
for( col in GE3_3P_logFits.WCU[,-1]){
  temp =  FoldCVGP(col,
                   myGD20_pruneCUfilt[which(missingGE3Plog),-1]-1, #convert to correct format
                   datasplit = .8,
                   trait_name =paste0('GELogistic_', colnames(GE3_3P_logFits.WCU)[counter],'_withVND'))
  temp2 =  FoldCVGP(col[which(GE3_3P_logFits.WCU$taxa %in% Mkk3taxa$taxa)],
                    (myGD20_pruneCUfilt[which(missingGE3Plog),] %>%
                       filter(taxa %in% Mkk3taxa$taxa))[,-1]-1, #convert to correct format
                    datasplit = .8,
                    trait_name =paste0('GELogistic_', colnames(GE3_3P_logFits.WCU)[counter],'_noVND'))
  temp3 = FoldCVGPSamePopSize(col,
                              myGD20_pruneCUfilt[which(missingGE3Plog),-1]-1, #convert to correct format
                              datasplit = .8,
                              trait_name =paste0('GELogistic_', colnames(GE3_3P_logFits.WCU)[counter],'_withVNDSamePop'),
                              sample_size = length(col[which(GE3_3P_logFits.WCU$taxa %in% Mkk3taxa$taxa)]))
  
  GPGElogfits =rbind(GPGElogfits,temp,temp2, temp3)
  
  counter = counter+1
}
# GP GE3 FPCA #####
GPGEFPCA =data.frame()
GE3.FPCA.WCU = GE31920FPCA.ouputs.longer %>% filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6'))%>%
  mutate(name = mapvalues(name, from = c('TimeTo90','TimeTo95'), to = c('fTimeTo90', 'fTimeTo95')))%>%
  filter(name != 'FPC4') %>%
  pivot_wider(names_from = 'name', values_from = 'value',) 
missingGE3FPCA= myGD20_pruneCUfilt$taxa %in% GE3.FPCA.WCU$taxa
counter = 2
for( col in GE3.FPCA.WCU[,-1]){
  temp =  FoldCVGP(col,
                   myGD20_pruneCUfilt[which(missingGE3FPCA),-1]-1, #convert to correct format
                   datasplit = .8,
                   trait_name =paste0('GEFPCA_',colnames(GE3.FPCA.WCU)[counter],'_withVND'))
  temp2 =  FoldCVGP(col[which(GE3.FPCA.WCU$taxa %in% Mkk3taxa$taxa)],
                    (myGD20_pruneCUfilt[which(missingGE3FPCA),] %>%
                       filter(taxa %in% Mkk3taxa$taxa))[,-1]-1, #convert to correct format
                    datasplit = .8,
                    trait_name =paste0('GEFPCA_',colnames(GE3.FPCA.WCU)[counter],'_noVND'))
  temp3 = FoldCVGPSamePopSize(col,
                              myGD20_pruneCUfilt[which(missingGE3FPCA),-1]-1, #convert to correct format
                              datasplit = .8,
                              trait_name = paste0('GEFPCA_',colnames(GE3.FPCA.WCU)[counter],'_withVNDSamePop'),
                              sample_size = length(col[which(GE3.FPCA.WCU$taxa %in% Mkk3taxa$taxa)]))
  
  GPGEFPCA =rbind(GPGEFPCA,temp,temp2,temp3)
  counter = counter+1
}

# GP GI3 Logistic ####
GPGIlogfits = data.frame()
GI3_4P_logFits.WCU = GI3_4P_logFits.w %>%
  filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6'))%>%
  mutate(TimeTo5.0 = as.numeric(ifelse(TimeTo5.0<250,TimeTo5.0,NA))) %>%
  select(!TimeTo5.6)
missingGI4Plog = myGD20_pruneCUfilt$taxa %in% GI3_4P_logFits.WCU$taxa
counter = 2
for( col in GI3_4P_logFits.WCU[,-1]){
  temp =  FoldCVGP(col,
                   myGD20_pruneCUfilt[which(missingGI4Plog),-1]-1, #convert to correct format
                   datasplit = .8,
                   trait_name =paste0('GILogistic_', colnames(GI3_4P_logFits.WCU)[counter],'_withVND'))
  temp2 =  FoldCVGP(col[which(GI3_4P_logFits.WCU$taxa %in% Mkk3taxa$taxa)],
                    (myGD20_pruneCUfilt[which(missingGI4Plog),] %>%
                       filter(taxa %in% Mkk3taxa$taxa))[,-1]-1, #convert to correct format
                    datasplit = .8,
                    trait_name =paste0('GILogistic_', colnames(GI3_4P_logFits.WCU)[counter],'_noVND'))
  temp3 =  FoldCVGPSamePopSize(col,
                               myGD20_pruneCUfilt[which(missingGI4Plog),-1]-1, #convert to correct format
                               datasplit = .8,
                               trait_name = paste0('GILogistic_', colnames(GI3_4P_logFits.WCU)[counter],'_withVNDSamePop'),
                               sample_size = length(col[which(GI3_4P_logFits.WCU$taxa %in% Mkk3taxa$taxa)]))
  
  GPGIlogfits =rbind(GPGIlogfits,temp,temp2, temp3)
  
  counter = counter+1
}

# GP GI3 FPCA ####
GPGIFPCA = data.frame()
GI3.FPCA.WCU = GI31920FPCA.ouputs.longer %>% filter(taxa %in% CULines$taxa &taxa %nin% c('P_5_3','P_8_6')) %>%
  mutate(name = mapvalues(name, from = c('TimeTo5.0','TimeTo5.6'),
                          to = c('fTimeTo5.0','fTimeTo5.6'))) %>% 
  pivot_wider(values_from = value, names_from = name) %>% select(!FPC4)
missingGI3FPCA= myGD20_pruneCUfilt$taxa %in% GI3.FPCA.WCU$taxa
counter = 2
for( col in GI3.FPCA.WCU[,-1]){
  temp =  FoldCVGP(col,
                   myGD20_pruneCUfilt[which(missingGI3FPCA),-1]-1, #convert to correct format
                   datasplit = .8,
                   trait_name =paste0('GIFPCA_', colnames(GI3.FPCA.WCU)[counter],'_withVND'))
  temp2 =  FoldCVGP(col[which(GI3.FPCA.WCU$taxa %in% Mkk3taxa$taxa)],
                    (myGD20_pruneCUfilt[which(missingGI3FPCA),] %>%
                       filter(taxa %in% Mkk3taxa$taxa))[,-1]-1, #convert to correct format
                    datasplit = .8,
                    trait_name =paste0('GIFPCA_', colnames(GI3.FPCA.WCU)[counter],'_noVND'))
  temp3 = FoldCVGPSamePopSize(col,
                              myGD20_pruneCUfilt[which(missingGI3FPCA),-1]-1, #convert to correct format
                              datasplit = .8,
                              trait_name = paste0('GIFPCA_', colnames(GI3.FPCA.WCU)[counter],'_withVNDSamePop'),
                              sample_size = length(col[which(GI3.FPCA.WCU$taxa %in% Mkk3taxa$taxa)]))
  
  GPGIFPCA =rbind(GPGIFPCA,temp, temp2, temp3)
  counter = counter+1
}

################################
# setwd('GWAS_OUTPUTS/')
# load("GE31920PerTP.GWAS.S.RData")
# load("GE3FPCA1920.GWAS.S.RData")
# load("GE3LogisticFits1920.GWAS.S.RData")
# load("GI31920PerTP.GWAS.S.RData")
# load("GI3FPCA1920.GWAS.S.RData")
# load("GI3LogisticFits1920.GWAS.S.RData")
# setwd(rprojroot::find_rstudio_root_file())

### PLOTS ######
### Figure 1 GE Recovered curve plots for FPCA and Logistic, and parameters histograms #########
GElogPcounts = GE3_3P_logFits.W %>% 
  pivot_longer(cols = 2:8) %>% filter(!is.na(value))%>% group_by(name) %>% 
  summarise(n = n()) %>% mutate(FacetLabel = paste0(name,': n=', n)) %>%
  merge(.,GE3_3P_logFits.W %>% pivot_longer(cols = 2:8),by = 'name') %>%
  mutate(ModelType = 'Logistic',
         ParamFiltered = taxa %in% (GE3_3P_logFits.W %>% filter(Lower>0.8))$taxa)
GEfpcaPcounts = GE31920FPCA.ouputs.longer %>% group_by(name) %>% filter(!(is.na(value)))%>%summarise(n = n()) %>% 
  mutate(FacetLabel = paste0(name,': n=', n)) %>% merge(.,GE31920FPCA.ouputs.longer, by ='name') %>% mutate(ModelType = 'FPCA')

GEfpcaRecoveredCurves = GE31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GE31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GE3', names_to = 'taxa') %>% select(taxa,time,GE3) %>% 
  mutate(ModelType = 'FPCA fits') %>% 
  ggplot(aes(x = time, y = GE3, group = taxa))+geom_line(color ='darkgrey')+
  facet_grid(cols = vars(ModelType))+
  geom_point(data = All.BluesGE31920%>% filter(taxa %nin% GE3logTaxatoFilter$taxa), color ='blue', size =.9)+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GE')

GElogisticRecoveredCurves = GE3fittedCurves %>% 
  filter(taxa %nin% GE3logTaxatoFilter$taxa) %>% 
  mutate(ModelType = 'Logistic fits',
         ParamFiltered = taxa %in% (GE3_3P_logFits.W %>% filter(Lower>0.8))$taxa) %>%
  ggplot(aes(x = time, y = GE3estimate, group = taxa))+geom_line(color = 'darkgrey')+
  facet_grid(cols = vars(ModelType)) +
  geom_point(data = All.BluesGE31920%>% filter(taxa %nin% GE3logTaxatoFilter$taxa),inherit.aes = F, 
             aes(x = time, y = GE3), color ='blue', size =.9)+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GE')

GECurveandParamPlot = GEfpcaRecoveredCurves+ 
  ggplot(GEfpcaPcounts %>%
           mutate(FacetLabel= mapvalues(FacetLabel, from = c('TimeTo90: n=484','TimeTo95: n=484'),to = c('fTimeTo90: n=484','fTimeTo95: n=484'))),
         aes(x = value))+geom_histogram()+facet_wrap(vars(FacetLabel),scales = 'free')+theme_bw(base_size = 8)+
  labs(subtitle = 'FPCA FPCs and derived values')+
  GElogisticRecoveredCurves+
  ggplot(GElogPcounts,aes(x = value))+geom_histogram()+
  facet_wrap(vars(FacetLabel), scales = 'free')+
  theme_bw(base_size = 8)+theme(legend.position = 'none')+
  labs(subtitle = 'Logistic parameter and derived values')+
  plot_layout(design = c('AABBB \n CCDDD'), heights = c(2,3))+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 8, hjust = 0, vjust = 0))
GECurveandParamPlot
setwd(rprojroot::find_rstudio_root_file())

pdf('Figure1_GECurveAndParameters.pdf',width =8,height = 6)
GECurveandParamPlot
dev.off()

### figure 2 GE Functional Manhattan Plot #####

GEfunctionalModelingManhattan = rbind(GE3FPCATopMarkers,GE3LogisticTopMarkers,
                                      filtGE3LogisticTopMarkers) %>%  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait),
         ModelType = mapvalues(ModelType, from = c('GE3Logistic','GE3FPCA','FiltGE3Logistic'),
                               to=c('Logistic','FPCA','Filtered Logistic'))) %>%
  mutate(trait = factor(trait, levels = c('FPC1','FPC2','FPC3','TimeTo90','TimeTo95','Centering',
                                          'Lower','Rate','rTimeTo90','rTimeTo95'))) %>%
  ggplot(aes(x = ordinal, y = log10PVal, shape = ModelType))+
  geom_point(aes(color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAT1')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+
  geom_hline(yintercept = -log10(5e-5))+
  ylim(0,30)+
  theme_bw()+labs(shape = 'Time series\nmodel',color = 'Parameter') +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1))

pdf('Figure2_GEFunctionalGWAresults.pdf', 8,5.5)
GEfunctionalModelingManhattan
dev.off()


### Figure 3 GI recovered curves and parameters plotted for FPCA and logistic curves #######
GIlogPcounts = GI3_4P_logFits.w %>% select(!c(TimeTo5.6)) %>%
  mutate(TimeTo5.0 = as.numeric(ifelse(TimeTo5.0<250,TimeTo5.0,NA)))%>%
  pivot_longer(cols = 2:9) %>% group_by(name)%>%
  filter(!is.na(value)) %>% summarise(n = n())%>% mutate(FacetLabel = paste0(name,': n=', n)) %>% 
  merge(.,GI3_4P_logFits.w %>% select(!c(TimeTo5.6)) %>%
          mutate(TimeTo5.0 = as.numeric(ifelse(TimeTo5.0<250,TimeTo5.0,NA)))%>% pivot_longer(cols = 2:9),
        by = 'name')
GIfpcaPcounts = GI31920FPCA.ouputs.longer %>% group_by(name) %>% filter(!(is.na(value)))%>%summarise(n = n()) %>% 
  mutate(FacetLabel= mapvalues(paste0(name,': n=', n), from = c('TimeTo5.0: n=478','TimeTo5.6: n=398'),
                               to = c('fTimeTo5.0: n=478','fTimeTo5.6: n=398'))) %>% 
  merge(.,GI31920FPCA.ouputs.longer, by ='name') %>% 
  mutate(ModelType = 'FPCA')

GIfpcaRecoveredCurves = GI31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GI31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GI3', names_to = 'taxa') %>%
  select(taxa,time,GI3) %>% mutate(ModelType = 'FPCA fits') %>% 
  ggplot(aes(x = time, y = GI3, group = taxa))+geom_line(color ='darkgrey')+
  facet_grid(cols = vars(ModelType))+
  geom_point(data = all1920GIBlues %>% filter(taxa %nin% GI3logTaxatoFilter$taxa), color ='blue', size = 0.9)+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GI')

GIlogisticRecoveredCurves = GI3fittedCurves %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% 
  mutate(ModelType = 'Logistic fits') %>%
  ggplot(aes(x = time, y = GI3estimate, group = taxa))+geom_line(color ='darkgrey')+
  facet_grid(cols = vars(ModelType)) +
  geom_point(data = all1920GIBlues%>% filter(taxa %nin% GI3logTaxatoFilter$taxa),
             inherit.aes = F, aes(x = time, y = GI3), color ='blue', size = 0.9)+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GI')

GICurveandParamPlot = 
  GIfpcaRecoveredCurves+ 
  ggplot(GIfpcaPcounts,aes(x = value))+geom_histogram()+facet_wrap(vars(FacetLabel),scales = 'free')+theme_bw(base_size = 8)+
  labs(subtitle = 'FPCA FPCs and derived values')+
  GIlogisticRecoveredCurves+
  ggplot(GIlogPcounts,aes(x = value))+geom_histogram()+facet_wrap(vars(FacetLabel), scales = 'free')+theme_bw(base_size = 8)+
  labs(subtitle = 'Logistic parameter and derived values')+
  plot_layout(design = c('AABBB \n CCDDD'), heights = c(2,3))+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 8, hjust = 0, vjust = 0))

GICurveandParamPlot
pdf('Figure3_GICurvesAndParameters.pdf', width = 8,height = 6)
GICurveandParamPlot
dev.off()

### Figure 4 GI Functional Manhattan #####
GIfunctionalModelingManhattan = rbind(GI3FPCATopMarkers,
                                      GI3LogisticTopMarkers,
                                      filtGI3LogisticTopMarkers) %>%  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait),
         ModelType = mapvalues(ModelType, from = c('GI3Logistic','GI3FPCA','FiltGI3Logistic' ), 
                               to=c('Logistic','FPCA', 'Filtered Logistic'))) %>%
  mutate(trait = factor(trait, levels = c('FPC1','FPC2','FPC3','TimeTo5.0','TimeTo5.6','Centering',
                                          'Lower','Rate','Upper','DeltaGI','DeltaGI90','DeltaGI95'))) %>%
  ggplot(aes(x = ordinal, y = log10PVal, shape = ModelType))+
  geom_point(aes(color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAT1')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+
  geom_hline(yintercept = -log10(5e-5))+
  ylim(0,30)+
  theme_bw()+labs(shape = 'Time series\nmodel',color = 'Parameter') +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1))

pdf('Figure4_GIFunctionalGWAresults.pdf', 8.00,6.00)
GIfunctionalModelingManhattan
dev.off()

### Figure 5 Results of GP on all parameters based on training sets #######
pdf('Figure5_SimpleGP.pdf', 8.00,6.00)
rbind(GPGEFPCA %>% mutate(trait = gsub(pattern = 'GE', replacement = '',x =trait)),
      GPGElogfits%>% mutate(trait = gsub(pattern = 'GE', replacement = '',x =trait)),
      GPPerTP %>% filter(substr(trait,10,12)=='GE3')%>% 
        mutate(trait = gsub(pattern = '3_', replacement = '_',x = trait))) %>% 
  group_by(trait, fold) %>% 
  summarize('Prediction accuracy' = mean(correlation)) %>% ungroup(trait, fold) %>% 
  separate(trait, sep = '_', into =c('xFacet', 'xMarker', 'MKK3included')) %>% mutate(yfacetLabel = 'GE') %>%
  mutate(xMarkerf = factor(xMarker,levels = c('FPC1','FPC2', 'FPC3','fTimeTo90','fTimeTo95','Centering','Lower','Rate',
                                              'TimeTo90', 'TimeTo95','rTimeTo90','rTimeTo95','TP1GE','TP2GE','TP3GE','TP4GE','TP5GE','TP6GE','TP7GE' )),
         Model = 'Model framework',
         xFacet = mapvalues(xFacet, from ='perTP', to = 'per Time Point'),
         MKK3included = factor(MKK3included, levels = c('withVND','noVND', 'withVNDSamePop'))) %>%
  ggplot(aes(x = `Prediction accuracy`, y = xMarkerf,  fill = MKK3included))+
  geom_boxplot()+
  facet_nested(Model+xFacet~yfacetLabel, scales = 'free',space = 'free_y')+
  theme_bw()+ylab(NULL)+scale_y_discrete(limits=rev)+
  xlim(-.2,1) + geom_vline(xintercept = 0)+
  scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb'),name ='MKK3 allele\nexcluded',
                    labels = c('None-full',expression(MKK3[N^{"*"}]), 'None-STPS'))+
  theme(legend.text.align = 0)+
  
  rbind(GPGIlogfits%>% mutate(trait = gsub(pattern = 'GI', replacement = '',x =trait)),
        GPGIFPCA %>% mutate(trait = gsub(pattern = 'GI', replacement = '',x =trait)),
        GPPerTP%>%filter(substr(trait,10,12)=='GI3') %>% 
          mutate(trait = gsub(pattern = '3_', replacement = '_',x = trait)),
        GPPHSCU) %>%
  group_by(trait, fold) %>% 
  summarize('Prediction accuracy' = mean(correlation)) %>% ungroup(trait, fold) %>% 
  separate(trait, sep = '_', into =c('xFacet', 'xMarker', 'MKK3included')) %>% 
  mutate(xMarker = mapvalues(xMarker, from =c('Delta','Delta90','Delta95'),
                             to = c('DeltaGI', 'DeltaGI90','DeltaGI95')),
         yfacetLabel = 'GI',
         xMarkerf = factor(xMarker,
                           levels = c('FPC1','FPC2', 'FPC3','fTimeTo5.0','fTimeTo5.6','Centering','Lower','Rate', 'Upper',
                                      'TimeTo5.0','DeltaGI90','DeltaGI95','DeltaGI','TP1GI','TP2GI','TP3GI','TP4GI',
                                      'TP5GI','TP6GI','TP7GI','PHS' )),
         Model = 'Model framework',
         xFacet = mapvalues(xFacet, from ='perTP', to = 'per Time Point'),
         MKK3included = factor(MKK3included, levels = c('withVND', 'noVND','withVNDSamePop'))) %>%
  arrange(MKK3included) %>%
  ggplot(aes(x = `Prediction accuracy`, y = xMarkerf, fill = MKK3included))+
  geom_boxplot()+
  facet_nested(Model+xFacet~yfacetLabel, scales = 'free',space = 'free_y')+
  theme_bw()+ ylab(NULL)+scale_y_discrete(limits=rev)+
  xlim(-.2,1) + geom_vline(xintercept = 0)+
  scale_fill_manual(values =  c('#66c2a5','#fc8d62','#8da0cb'),name ='MKK3 allele\nexcluded',
                    labels = c('None-full',expression(MKK3[N^{"*"}]), 'None-STPS'))+
  theme(legend.text.align = 0)+
  plot_layout(ncol = 2, guides = 'collect')+plot_annotation(tag_levels = 'a')
dev.off()

### Figure 6 Leave one TP out cross validation ######
setwd(rprojroot::find_rstudio_root_file())
load("Output/FPCA_GEGI_predictions.RData")
load("Output/Logistic_GEGI_predictions.RData")

rectangles = FPCA_GEGI_predictions %>% select(TP, Trait, Dropped) %>% unique()%>%
  mutate(Top = 5, Bottom = -5, TraitLabel = paste0('Trait: ', Trait), TP = Dropped) %>%
  filter(Dropped !='None') %>% mutate(alpha = 0.05)

FPCALogisticPrediction = Logistic_GEGI_predictions %>%
  mutate(TraitLabel = paste0('Model: Logistic; Trait: ', Trait)) %>%
  rbind(FPCA_GEGI_predictions %>%
          mutate(TraitLabel = paste0('Model: FPCA; Trait: ', Trait))) %>% 
  mutate(TraitLabel = factor(TraitLabel,levels = c('Model: FPCA; Trait: GE',
                                                   'Model: Logistic; Trait: GE',
                                                   'Model: FPCA; Trait: GI',
                                                   'Model: Logistic; Trait: GI')),
         MethodLabel = factor(mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                        to = c('Training Set\nAccuracy',
                                               'Testing Set\nGPP Accuracy', 'Testing Set\nGPT Accuracy')),
                              levels = c('Training Set\nAccuracy',
                                         'Testing Set\nGPT Accuracy','Testing Set\nGPP Accuracy'),
                              ordered = TRUE),
         TPModel = paste0(TP,':',Model))
Rectangles2 = rbind(rectangles %>%
                      mutate(TraitLabel = paste0('Model: Logistic; Trait: ', Trait),
                             Model = 'Logistic'),
                    rectangles %>%
                      mutate(TraitLabel = paste0('Model: FPCA; Trait: ', Trait),
                             Model = 'FPCA')) %>%
  mutate(TraitLabel = factor(TraitLabel,levels = c('Model: FPCA; Trait: GE',
                                                   'Model: Logistic; Trait: GE',
                                                   'Model: FPCA; Trait: GI',
                                                   'Model: Logistic; Trait: GI')),
         TPModel = paste0(TP,':',Model))

library(ggh4x)
pdf('Figure6_LeaveOneTPoutGP.pdf', 10,8)
FPCALogisticPrediction %>% mutate('TP Dropped' = 'Time point masked') %>%
  ggplot(aes(x = TP, y = correlation, color = MethodLabel))+
  geom_boxplot(outlier.size = 0.5) +
  facet_nested(`TP Dropped`+Dropped  ~ TraitLabel)+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a'))+
  geom_hline(yintercept = 0, color = 'black') +
  xlab('Time point of prediction')+
  labs(color = 'Prediction \nMethod') +
  geom_tile(data = Rectangles2 %>% mutate('TP Dropped' = 'Time point masked'),
            inherit.aes = FALSE,
            aes(x = TP, fill = NA,y = 0), alpha = 0.05,
            height = 2, width = 1, fill = 'lightgrey') +
  theme_bw(base_size = 10) +
  ylab('Prediction accuracy')
dev.off()


### Figure 7 Time To thresholds vs PHS and other things ######
setwd(rprojroot::find_rstudio_root_file())
load("PhenotypeData/ProcessedData/PHS_BLUEs_GGS1920.RData")
PopulationKey = read.csv(file = 'GenotypeData/PopulationKey.csv') %>% select(!X)
CULines =PopulationKey %>% filter(Population %in% c("check","C1G","C2G","C1P","base"))

#PHS with ftimeto5.0 and ftimeto95 in ONLY THE CU LINES
PHS.blues %>% merge(.,
                    all_BLUE %>% select(TP, taxa, GE3,GI3scale)%>%filter(TP =='TP1'),
                    by = 'taxa')%>%
  rename('GE TP1' =GE3, 'GI TP1' = GI3scale) %>%
  merge(.,VectorofTimeto95[c('taxa',"TimeTo95")], by = 'taxa') %>%
  merge(.,VectorofTimeto5.0[c('taxa','TimeTo5.0')], by = 'taxa') %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342', 'AlaAT_L214F')], by = 'taxa',all.x = TRUE) %>% 
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>% 
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         AlaAT = round(AlaAT_L214F),
         AlaATcode = mapvalues(AlaAT,from =c(0,2,1), to = c('N','D','Het')),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), to= c('Dormant','Non-Dormant','Very \n Non-Dormant'))) %>% 
  rename(fTimeTo95 = TimeTo95, fTimeTo5.0 = TimeTo5.0) %>%
  filter(taxa %in% CULines$taxa) %>%
  gatherpairs(PHS, 'GE TP1', 'GI TP1',fTimeTo95,fTimeTo5.0) %>%
  ggplot(aes(x = .xvalue, y = .yvalue, color = SD2Code))+
  geom_point()+
  scale_colour_manual(values = c('#66c2a5','#fc8d62','#8da0cb'),name = "*HvMKK3* \nallele", 
                      labels = expression(MKK3[D],MKK3[N],MKK3[N^{"*"}]))+theme_bw() +
  theme(legend.title = element_markdown())+
  facet_grid(cols = vars(.xkey), rows = vars(.ykey), scales = 'free',switch = "y")+
  xlab('Values')+ylab('Values')

#lets try with ggally to get what I want...
#This is mostly what I want expect for the markdown problems with the labeller. 
pdf('Figure7_PHScorrelations.pdf', 10.00,6.00)
PHS.blues %>% merge(.,
                    all_BLUE %>% select(TP, taxa, GE3,GI3scale)%>%filter(TP =='TP1'),
                    by = 'taxa')%>%
  rename('GETP1' =GE3, 'GITP1' = GI3scale) %>%
  merge(.,VectorofTimeto95[c('taxa',"TimeTo95")], by = 'taxa') %>%
  merge(.,VectorofTimeto5.0[c('taxa','TimeTo5.0')], by = 'taxa') %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342', 'AlaAT_L214F')], by = 'taxa',all.x = TRUE) %>%
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>%
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), 
                            to =c('MKK3 D','MKK3 N','MKK3 N*'))) %>%
  rename(fTimeTo95 = TimeTo95, fTimeTo5.0 = TimeTo5.0) %>%
  filter(taxa %in% CULines$taxa) %>%
  ggpairs(data = ., 
          columns = c('PHS', 'GETP1', 'GITP1','fTimeTo95','fTimeTo5.0'),
          ggplot2::aes(color = SD2Code))+theme_bw()+theme(axis.text = element_text(size = 6))
dev.off()

### Sup Table 1 and 2 significant per TP and Functional Markers for the tables ######
# per TP
rbind(GI3perTPTopMarkers,GE3perTPTopMarkers) %>%  filter(maf > 0.07 &P.value<5e-5) %>% 
  mutate(log10PVal = round(log10PVal,2)) %>%
  select(ModelType,trait, SNP,Chromosome, Position,maf, log10PVal) 

#for functional analysis
SigMarkers = rbind(GI3FPCATopMarkers, 
                   GE3FPCATopMarkers,
                   GI3LogisticTopMarkers,
                   GE3LogisticTopMarkers,
                   filtGE3LogisticTopMarkers,
                   filtGI3LogisticTopMarkers) %>%  filter(maf >0.07 & P.value < 5e-5) %>%
  arrange(Chromosome, Position)
View(SigMarkers %>% select(ModelType,trait,SNP, Chromosome, Position,maf, P.value))
SigMarkers %>%
  arrange(ModelType, trait, Chromosome, Position, P.value) %>%
  mutate(log10PVal = round(log10PVal,2)) %>%
  select(ModelType,trait, SNP,Chromosome, Position, maf, log10PVal)%>%
  merge(., rbind(
    GE3_3P_logFits.W %>% 
      pivot_longer(cols = 2:8) %>% filter(!is.na(value))%>% group_by(name) %>% 
      summarise(n = n()) %>%  mutate(ModelType = 'GE3Logistic') %>% rename(trait = name),
    GE3_3P_logFits.W %>% filter(Lower<0.8) %>%
      pivot_longer(cols = 2:8) %>% filter(!is.na(value))%>% group_by(name) %>% 
      summarise(n = n()) %>%  mutate(ModelType = 'FiltGE3Logistic') %>% rename(trait = name),
    GE31920FPCA.ouputs.longer %>% group_by(name) %>% filter(!(is.na(value)))%>%summarise(n = n()) %>% 
      mutate(ModelType = 'GE3FPCA') %>% rename(trait = name),
    GI3_4P_logFits.w %>% 
      mutate(TimeTo5.0 = as.numeric(ifelse(TimeTo5.0<250,TimeTo5.0,NA)))%>%
      pivot_longer(cols = 2:10) %>% group_by(name)%>%
      filter(!is.na(value)) %>% summarise(n = n())%>%mutate(ModelType = 'GI3Logistic') %>% 
      rename(trait = name),
    GI3_4P_logFits.w %>% filter(Lower<5.0) %>%  
      mutate(TimeTo5.0 = as.numeric(ifelse(TimeTo5.0<250,TimeTo5.0,NA)))%>%
      pivot_longer(cols = 2:10) %>% group_by(name)%>%
      filter(!is.na(value)) %>% summarise(n = n())%>% mutate(ModelType = 'FiltGI3Logistic') %>% 
      rename(trait = name),
    GI31920FPCA.ouputs.longer %>% group_by(name) %>% filter(!(is.na(value)))%>%summarise(n = n()) %>% 
      rename(trait = name) %>% mutate(ModelType = 'GI3FPCA')), by = c('ModelType','trait'),all.X =TRUE)%>%
  mutate(ModelType = mapvalues(ModelType, 
                               from=c("GI3Logistic","GI3FPCA","GE3Logistic","FiltGE3Logistic","FiltGI3Logistic","GE3FPCA"),
                               to=c("GILogistic","GIFPCA","GELogistic","GELogisticFilt","GILogisticFilt","GEFPCA")))%>%
  rename(Chr = Chromosome, Trait = trait, '-log(p-value)' = log10PVal) %>%
  mutate('Marker Region' = ifelse(Chr==5 & Position>585246000 & Position < 600000000, 'SD2',
                                  ifelse(Chr==5 &Position >442000000 & Position <443000000,'SD1',NA)),
         'Gene Candidate' = mapvalues(`Marker Region`, from = c('SD1','SD2'), to = c('HvAlaAT1','HvMKK3'))) %>%
  tail(.,17)

uniqueMarkers = SigMarkers %>% select(SNP) %>% mutate(SNP = as.character(SNP)) %>% unique()
ModelList = c(rep("GIFPCA",5),
              rep("GEFPCA",5),
              rep("GILogistic",8),
              rep("GELogistic",7),
              rep("GELogisticFilt",7),
              rep("GILogisticFilt",8))
TraitList = c(GI3FPCA1920.GWAS.S.traits,GE3FPCA1920.GWAS.S.traits,GI3LogisticFits1920.GWAS.S.Traits,
              GE3LogisticFits1920.GWAS.S.Traits,filtGE3LogisticFits1920.GWAS.S.Traits,filtGI3LogisticFits1920.GWAS.S.Traits)
AllPvals = data.frame()
Countout = 1
counter = 1
for (i in c(GI3FPCA1920.GWAS.S,GE3FPCA1920.GWAS.S,GI3LogisticFits1920.GWAS.S,GE3LogisticFits1920.GWAS.S,
            FiltGE3LogisticFits1920.GWAS.S,FiltGI3LogisticFits1920.GWAS.S)){
  print(head(i))
  AllPvals = rbind(i %>% mutate(SNP = as.character(SNP))%>% filter(SNP %in% uniqueMarkers$SNP) %>% 
                     select(SNP,Chromosome, Position,maf, log10PVal) %>%
                     mutate(ModelType =ModelList[counter],
                            Trait = TraitList[counter])
                   ,AllPvals)
  counter =counter+1  
}
AllPvals %>% colnames()
AllPvals %>%select(!maf)%>%pivot_wider(names_from = c(ModelType,Trait), names_sep = '_', values_from = log10PVal)


### Sup figure 1 GE FPCA models FPC effects coloring ####
jpeg(filename = 'GEFPCscolored.jpg', 700,350, res = 120)
GE31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GE31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GE3', names_to = 'taxa') %>%
  select(taxa,time,GE3) %>% mutate(ModelType = 'FPC1') %>% 
  merge(GE31920.FPCA$PCs_withTaxa, all.x = TRUE, by = 'taxa') %>%
  ggplot(aes(x = time, y = GE3, group = taxa))+geom_line(aes(color = FPC1))+
  facet_grid(cols = vars(ModelType))+
  # labs(subtitle = 'FPCA fits')+
  # geom_point(data = All.BluesGE31920%>% filter(taxa %nin% GE3logTaxatoFilter$taxa), color ='blue', size =.9)+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GE') +
  GE31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GE31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GE3', names_to = 'taxa') %>%
  select(taxa,time,GE3) %>% mutate(ModelType = 'FPC2') %>% 
  merge(GE31920.FPCA$PCs_withTaxa, all.x = TRUE, by = 'taxa') %>%
  ggplot(aes(x = time, y = GE3, group = taxa))+geom_line(aes(color = FPC2))+
  facet_grid(cols = vars(ModelType))+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GE')
dev.off()

### Sup figure 2 GI FPCA models FPC effects coloring #####
jpeg(filename = 'GIFPCscolored.jpg', 700,350, res = 120)
GI31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GI31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GI3', names_to = 'taxa') %>%
  select(taxa,time,GI3) %>% mutate(ModelType = 'FPC1') %>% 
  merge(GI31920.FPCA$PCs_withTaxa, all.x = TRUE, by = 'taxa') %>%
  ggplot(aes(x = time, y = GI3, group = taxa))+geom_line(aes(color = FPC1))+
  facet_grid(cols = vars(ModelType))+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GI') +
  GI31920.FPCA$RecoveredCurves %>% as.data.frame() %>% mutate(time = GI31920.FPCA$phi.fun.df$time) %>%
  pivot_longer(cols = 1:484, values_to = 'GI3', names_to = 'taxa') %>%
  select(taxa,time,GI3) %>% mutate(ModelType = 'FPC2') %>% 
  merge(GI31920.FPCA$PCs_withTaxa, all.x = TRUE, by = 'taxa') %>%
  ggplot(aes(x = time, y = GI3, group = taxa))+geom_line(aes(color = FPC2))+
  facet_grid(cols = vars(ModelType))+
  theme_bw(base_size = 8) +xlab('Time')+ylab('GI')
dev.off()

### Sup Figure 3 per time point Manhattan plot #####
jpeg('perTPModeling.jpg', 800,400, res = 120)
rbind(GI3perTPTopMarkers,GE3perTPTopMarkers) %>%
  mutate(ModelType = mapvalues(ModelType, from = c('GE3 per Time Point','GI3 per Time Point'),
                               to =c('GE','GI') )) %>%
  filter(maf >0.07) %>%
  mutate(ModelTypeParam = paste0(ModelType,trait),
         trait = substr(trait,1,3)) %>%
  ggplot(aes(x = ordinal, y = log10PVal, shape = ModelType))+
  geom_point(aes(color = trait), size =2.5) +
  geom_vline(xintercept = ChromLines)+
  geom_vline(xintercept = c(9228,10771), color = 'red')+
  annotate(geom = 'text', x = 9228, y = 22, label = 'AlaAt')+
  annotate(geom = 'text', x = 10771, y = 18, label = 'MKK3')+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = breaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+
  geom_hline(yintercept = -log10(5e-5))+
  theme_bw()+labs(color = 'Time point',shape = 'Trait')+
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1))

dev.off()

### Sup Figure 4 GE Logistic Parameters correlations and colored by MKK3 status #########
jpeg(filename = 'GELogParamsScaterplot.jpg',700,500,res= 120 )
merge(GE3_3P_logFits.W%>%select(taxa, Lower,Rate, Centering) %>%
        pivot_longer(.,names_to = 'Xnames',values_to = 'Xvalues', cols = c(Lower,Rate,Centering)),
      GE3_3P_logFits.W%>%select(taxa, Lower,Rate, Centering), by = 'taxa')%>%
  pivot_longer(.,names_to = 'Ynames',values_to = 'Yvalues', cols = c(Lower,Rate,Centering)) %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342')], by = 'taxa',all.x = TRUE) %>% 
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>%
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), 
                            to = c('D','ND','VND'))) %>%
  ggplot(aes(x = Xvalues, y = Yvalues, color = SD2Code))+geom_point()+
  facet_grid(rows = vars(Ynames), cols = vars(Xnames), scales = 'free')+theme_bw() +
  scale_colour_manual(values = c('#66c2a5','#fc8d62','#8da0cb'),name = "*HvMKK3* \nallele", 
                      labels = expression(MKK3[D],MKK3[N],MKK3[N^{"*"}]))+
  theme(legend.title = element_markdown())+
  xlab('Values')+ylab('Values')
dev.off()


# set dplyr functions
select <- dplyr::select; rename <- dplyr::rename; mutate <- dplyr::mutate; 
summarize <- dplyr::summarize; arrange <- dplyr::arrange; slice <- dplyr::slice; filter <- dplyr::filter; recode<-dplyr::recode

# remove obs for setosa
setwd(rprojroot::find_rstudio_root_file())
jpeg(filename = 'GELogParamsGGallyPlot.jpg',700,500,res= 120 )
GE3_3P_logFits.W%>%select(taxa, Lower,Rate, Centering) %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342')], by = 'taxa',all.x = TRUE) %>% 
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>%
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), 
                            to = c('MKK3 D','MKK3 N','MKK3 N*'))) %>%
  ggpairs(data = ., columns = c('Lower','Rate','Centering'),
          ggplot2::aes(color = SD2Code))+theme_bw() 
# scale_color_manual(values=c('MKK3<sub>D</sub>'="#66c2a5",'MKK3<sub>N</sub>'='#fc8d62','MKK3<sub>N*</sub>'='#8da0cb',"Overall Corr"="black"))
dev.off()


### Sup Figure 5 Time to 95 and 90 GE logistic models png histograms KS tests  ######
testing = data.frame(taxa = myGD20_prune$taxa,
                     Chr6mk = paste0(round(myGD20_prune$`JHI-Hv50k-2016-408912`,0)),
                     SDHaplo = paste0(round(myGD20_prune$AlaAT_L214F,0),
                                      round(myGD20_prune$`JHI-Hv50k-2016-367342`,0),
                                      myGD20_prune$MKK3_E165Q)) %>%
  merge(VectorofTimeto95[c('taxa','TimeTo95')], by = 'taxa',all = TRUE) %>% rename(fTimeTo95 = TimeTo95) %>%
  merge(GE3_3P_logFits.W[c('taxa','TimeTo95')],by ='taxa', all = TRUE) %>%
  merge(VectorofTimeto90[c('taxa','TimeTo90')], by = 'taxa',all = TRUE) %>% rename(fTimeTo90 = TimeTo90) %>%
  merge(GE3_3P_logFits.W[c('taxa','TimeTo90')],by ='taxa', all = TRUE)
cor(as.numeric(testing$Chr6mk), testing$TimeTo90, use = 'complete.obs')
testing %>% select(TimeTo90, fTimeTo90, TimeTo95, fTimeTo95) %>% 
  pairs(.,lower.panel = upper.panel,upper.panel = panel.cor)

timetoDormbreaksum = testing %>% 
  pivot_longer(cols =c(TimeTo90, fTimeTo90, TimeTo95, fTimeTo95)) %>%
  group_by(name) %>%  filter(!is.na(value)) %>%
  summarize(Mean = mean(value),
            '  0th Percentile' = quantile(value)[1],
            ' 25th Percentile' = quantile(value)[2],
            ' 50th Percentile, \n Median' = quantile(value)[3],
            ' 75th Percentile' = quantile(value)[4],
            '100th Percentile' = range(value)[2]) %>% 
  pivot_longer(cols = !name, names_to = 'Statistic') %>%
  mutate(linetype = ifelse(Statistic == 'Mean','solid','dashed'))

jpeg('TimeToDormancyBreakHist.jpg', 700,500, res = 120)
testing %>%  pivot_longer(cols =c(TimeTo90, fTimeTo90, TimeTo95, fTimeTo95)) %>% 
  filter(value<100) %>%
  ggplot(aes(x= value)) +facet_wrap(vars(name))+geom_histogram()+
  geom_vline(data = timetoDormbreaksum,
             aes(xintercept = value, color = Statistic, linetype =linetype), size = 1.6) +
  xlim(0, 100) + guides(linetype = FALSE) + theme_bw() +xlab('Estimated time to threshold')
dev.off()

ks.test(x = testing$TimeTo90, y = testing$fTimeTo90, alternative = c('two.sided')) 
ks.test(x = testing$TimeTo95, y = testing$fTimeTo95, alternative = c('two.sided')) 
t.test(x = testing$TimeTo90, y = testing$fTimeTo90, alternative = c('two.sided')) 
t.test(x = testing$TimeTo95, y = testing$fTimeTo95, alternative = c('two.sided')) 



### Sup Figure 6 GI Logistic parameters Correlations and colored by MKK3 Status ######

jpeg(filename = 'GILogParamsScaterplot.jpg',700,500,res= 120 )
merge(GI3_4P_logFits.w%>%select(taxa, Lower,Rate, Centering, Upper) %>%
        pivot_longer(.,names_to = 'Xnames',values_to = 'Xvalues', cols = c(Lower,Rate,Centering,Upper)),
      GI3_4P_logFits.w%>%select(taxa, Lower,Rate, Centering, Upper), by = 'taxa')%>%
  pivot_longer(.,names_to = 'Ynames',values_to = 'Yvalues', cols = c(Lower,Rate,Centering,Upper)) %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342')], by = 'taxa',all.x = TRUE) %>% 
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>%
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), to= c('Dormant','Non-Dormant','Very \n Non-Dormant'))) %>%
  ggplot(aes(x = Xvalues, y = Yvalues, color = SD2Code))+geom_point()+
  facet_grid(rows = vars(Ynames), cols = vars(Xnames), scales = 'free')+theme_bw() +
  scale_colour_manual(values = c('#66c2a5','#fc8d62','#8da0cb'),name = "*HvMKK3* \nallele", 
                      labels = expression(MKK3[D],MKK3[N],MKK3[N^{"*"}]))+
  theme(legend.title = element_markdown())+
  xlab('Values')+ylab('Values')
dev.off()

jpeg(filename = 'GILogParamsGGallyPlot.jpg',700,500,res= 120 )
data.frame(taxa = GI3_4P_logFits.w$taxa,
           Rate = as.vector(GI3_4P_logFits.w$Rate),
           Upper = as.vector(GI3_4P_logFits.w$Upper),
           Centering = as.vector(GI3_4P_logFits.w$Centering),
           Lower = as.vector(GI3_4P_logFits.w$Lower)) %>%
  merge(., myGD20_prune[c('taxa','MKK3_E165Q','JHI-Hv50k-2016-367342')], by = 'taxa',all.x = TRUE) %>% 
  filter(MKK3_E165Q != 1 & `JHI-Hv50k-2016-367342` != 1)%>%
  mutate(SD2 = paste0(round(`JHI-Hv50k-2016-367342`),round(MKK3_E165Q)),
         MKK3_E165Q = paste(MKK3_E165Q),
         SD2Code= mapvalues(SD2, from = c('00','20','22'), 
                            to =c('MKK3 D','MKK3 N','MKK3 N*'))) %>%
  ggpairs(data = ., 
          columns = c('Lower','Rate','Centering','Upper'),
          ggplot2::aes(color = SD2Code))+theme_bw()
# scale_color_manual(values=c('MKK3<sub>D</sub>'="#66c2a5",'MKK3<sub>N</sub>'='#fc8d62','MKK3<sub>N/*</sub>'='#8da0cb',"Overall Corr"="black"))
dev.off()


### Sup Figure 7 Some examples of poor GI fits #####
jpeg('PoorGILogFits.jpg', 500,400, res = 120)
GI3fittedCurves %>% filter(taxa %nin% GI3logTaxatoFilter$taxa) %>% 
  mutate(ModelType = 'Logistic fits with rate < -10') %>% filter(taxa %in% (GI3_4P_logFits.w %>% filter(Rate < -10))$taxa)%>%
  ggplot(aes(x = time, y = GI3estimate, group = taxa, color = taxa))+geom_line(size = 1.5)+
  facet_grid(cols = vars(ModelType)) +
  geom_point(data = all1920GIBlues%>% filter(taxa %nin% GI3logTaxatoFilter$taxa)%>% 
               filter(taxa %in% (GI3_4P_logFits.w %>% filter(Rate < -10))$taxa),
             inherit.aes = F, aes(x = time, y = GI3,color = taxa), size = 1.5)+
  theme_bw(base_size = 10 ) +xlab('Time')+ylab('GI') +labs(color = 'Line')
dev.off()
### Distributions of fTimeto5.0 vs TimeTo5.0
jpeg('TimeTo5thresholdsfpcavslogsitic.jpg', 500,400, res = 120)
join(VectorofTimeto5.0 %>% select(taxa, TimeTo5.0) %>% rename(fTimeTo5.0 = TimeTo5.0),
     GI3_4P_logFits %>% filter(term == 'TimeTo5.0') %>% select(taxa, estimate) %>%
       filter(estimate<250) %>%
       rename(TimeTo5.0 = estimate)) %>%
  pivot_longer(cols = !taxa) %>%
  ggplot(aes(x = value))+geom_histogram()+
  facet_wrap(vars(name), ncol = 1)+
  geom_vline(data = join(VectorofTimeto5.0 %>% select(taxa, TimeTo5.0) %>% rename(fTimeTo5.0 = TimeTo5.0),
                         GI3_4P_logFits %>% filter(term == 'TimeTo5.0') %>% select(taxa, estimate) %>%
                           filter(estimate<250) %>%
                           rename(TimeTo5.0 = estimate)) %>%  
               pivot_longer(cols = !taxa)%>% filter(!is.na(value)) %>%
               group_by(name) %>% 
               summarize(Mean = mean(value),
                         '  0th Percentile' = quantile(value)[1],
                         ' 25th Percentile' = quantile(value)[2],
                         ' 50th Percentile;\nMedian' = quantile(value)[3],
                         ' 75th Percentile' = quantile(value)[4]) %>%
               pivot_longer(cols = !name, names_to = 'Statistic') %>%
               mutate(linetype = ifelse(Statistic == 'Mean','solid','dashed')),
             aes(xintercept = value, color = Statistic), size = 1.6)+
  xlab('Time to threshold') +theme_bw()
dev.off()




### LD between things.  #############
myGM20_prune %>% filter(Chromosome==6 & Position>472000000 & Position<474000000)
ld = myGD20_prune %>% select(AlaAT_L214F,MKK3_E165Q, `JHI-Hv50k-2016-367342`,`JHI-Hv50k-2016-408912`,
                             `JHI-Hv50k-2016-408918`,`JHI-Hv50k-2016-408820`,
                             `JHI-Hv50k-2016-408472`)
markerList = myGM20_prune %>% filter(SNP %in% c('AlaAT_L214F','MKK3_E165Q', 'JHI-Hv50k-2016-367342','JHI-Hv50k-2016-408912',
                                                'JHI-Hv50k-2016-408918','JHI-Hv50k-2016-408820',
                                                'JHI-Hv50k-2016-408472'))

df = ld

ld_heatmap=function(df){
  ld <- as.matrix(round(df,0))
  
  if(c(-1,3,4) %in% ld){
    ld[which(ld==3)]=2
    ld[which(ld==4)]=2
    ld[which(ld== -1)]=0
  }
  
  LD <- LD.Measures(donnees=ld,  na.presence=F)
  #LD$loc1=as.character(LD$loc1); LD$loc2=as.character(LD$loc2)
  r2 <- matrix(0, nrow=ncol(df), ncol=ncol(df))
  r2[lower.tri(r2, diag=FALSE)] <- LD$r2
  r2 <- t(r2)
  r2 <- as.data.frame(round(r2, 5))
  diag(r2) <- 1
  r2[lower.tri(r2)] = NA
  rownames(r2)=colnames(df); colnames(r2)=rownames(r2)
  r_2=melt(as.matrix(r2), na.rm=T)
  
  r_2  
  graphic = ggplot(r_2, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", name="r2") +
    theme_classic() +
    geom_text(aes(label = value))+theme_bw()#+    ggtitle(paste("LD r2 from",colnames(r2)[1],"-", colnames(r2)[length(colnames(r2))], sep=" " ))
  return(graphic)
}

ld_heatmap(ld) + labs(title = 'LD between SD1, SD2, and Chromosome 6H Markers')+
  xlab('')+ylab('') +theme(axis.text.x = element_blank())

### Random Correlations ###########

merge(VectorofTimeto5.0[c('taxa','TimeTo5.0')],
      VectorofTimeto5.6[c('taxa','TimeTo5.6')], by ='taxa', all.x = TRUE) %>%
  select(TimeTo5.0,TimeTo5.6) %>% cor(.,use = 'complete.obs')

cor(GI3_4P_logFits.w$DeltaGI90,GI3_4P_logFits.w$DeltaGI95, use = 'complete.obs')
merge(GE3_3P_logFits.W[c('taxa','TimeTo90')],
      VectorofTimeto90[c('taxa','TimeTo90')], by = 'taxa', all = TRUE) %>%
  select(TimeTo90.x,TimeTo90.y) %>% cor(use = 'complete.obs', method = 'spearman')
merge(GE3_3P_logFits.W[c('taxa','TimeTo90')],
      VectorofTimeto90[c('taxa','TimeTo90')], by = 'taxa', all = TRUE) %>%
  select(TimeTo90.x,TimeTo90.y) %>% plot()

merge(GE3_3P_logFits.W[c('taxa','TimeTo95')],
      VectorofTimeto95[c('taxa','TimeTo95')], by = 'taxa', all = TRUE) %>%
  select(TimeTo95.x,TimeTo95.y) %>% cor(use = 'complete.obs', method = 'spearman')

merge(GE3_3P_logFits.W[c('taxa','TimeTo95')],
      VectorofTimeto95[c('taxa','TimeTo95')], by = 'taxa', all = TRUE) %>%
  select(TimeTo95.x,TimeTo95.y) %>% plot()

cor(myGD20_prune$`JHI-Hv50k-2016-367342`, myGD20_prune$`JHI-Hv50k-2016-366325`)

GI3_4P_logFits.w %>% select(taxa,TimeTo5.0) %>% mutate(TimeTo5.0 = ifelse(TimeTo5.0 >250,NA,TimeTo5.0))%>%
  merge(VectorofTimeto5.0, by = 'taxa') %>%
  select(TimeTo5.0.x, TimeTo5.0.y) %>% cor(.,use = 'complete.obs')
GI3_4P_logFits.w %>% select(taxa,TimeTo5.0) %>% mutate(TimeTo5.0 = ifelse(TimeTo5.0 >250,NA,TimeTo5.0))%>%
  merge(VectorofTimeto5.0, by = 'taxa') %>%
  select(TimeTo5.0.x, TimeTo5.0.y) %>% plot()

ks.test(x = GI3_4P_logFits.w %>% select(TimeTo5.0) %>% mutate(TimeTo5.0 = ifelse(TimeTo5.0 >250,NA,TimeTo5.0)),
        y = VectorofTimeto5.0$TimeTo5.0, alternative = c('two.sided')) 




