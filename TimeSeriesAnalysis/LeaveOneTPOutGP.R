# Genomic prediction for the time point masked and full data sets 
# FunctionFitting_GWAS_SimpleGP.R MUST BE RUN BEFORE RUNNING THIS SCRIPT!!
# Caution: It takes a few hours to run!
library(Hmisc);library(drc);library(ggplot2);library(readxl);library(reshape2);
library(patchwork);library(rrBLUP);
library(fda) ; library(magic) ;
library(plyr);library(dplyr);
library(knitr);library(tidyr)

# FPCA Function for full and cases where a TP is dropped. #####
GermFPCAcomplexGP = function(GIlines = all1920GIBlues.FPCAInput %>%  filter(taxa %in% CULines$taxa & taxa %nin% c('P_5_3','P_8_6'))%>%
                             select(taxa,time,GI3), #I know is weird, but all the inputs for this need to be labels for GI regarless of if it is for GE
                           myGD = myGD20_pruneCUfilt,
                           numfolds = 50,
                           datasplit = .80,
                           Truth = all1920GIBlues %>% select(taxa, time,TP, GI3)){
  GIlines = GIlines %>% arrange(taxa, time)
  checkorder = sum(GIlines %>% filter(time ==0) %>% select(taxa) == 
                     as.character(((myGD %>% filter(taxa %in% GIlines$taxa) %>% arrange(taxa))$taxa)))
  myGD = (myGD %>% filter(taxa %in% GIlines$taxa) %>% arrange(taxa))[,-1]-1
  Results = data.frame()
  set.seed(1)
  num_entries = length(unique(GIlines$taxa))
  Training_size = round(datasplit*num_entries,0)
  numObs = as.integer((GIlines %>%group_by(taxa) %>% summarise(n()))[1,2])
  if(checkorder != num_entries){
    stop('myGD and the phenotypes are not aligning correctly, ensure order')
  }
  for (ii in 1:numfolds){
    start_time <- Sys.time()
    trainingnums = sort(sample(1:num_entries,Training_size))
    testingnums = setdiff(1:num_entries,trainingnums)

    trainingSet = GIlines %>% filter(taxa %in% unique(GIlines$taxa)[trainingnums])
    testingSet = GIlines %>% filter(taxa %in% unique(GIlines$taxa)[testingnums])
    marker_train = myGD[trainingnums,]
    marker_test = myGD[testingnums,]
    FPCA = FPCA_function(dfTaxaTraitTime = trainingSet,
                  Trait = 'GI3', #Trait name must be entered as a character ie Trait = 'GE3'
                  NumKnots = 1, # NumKnots is the number of interior knots to be fitted
                  order = 3, # Order is the dergree of the polynomial to be fit to the data.
                  NumObsevationsPerLine =numObs )  
    numPcs = sum(FPCA$v1>0)
    PredParams = data.frame(taxa = unique(GIlines$taxa)[testingnums])
    counter = 2
    for (col in FPCA$PCs_withTaxa[,c(2:(numPcs+1))]){
      trained.Model = mixed.solve(col, Z = marker_train, K = NULL, SE =FALSE)
      PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u) + as.numeric(trained.Model$beta)
      print(trained.Model$beta)
      PredParams[[(colnames(FPCA$PCs_withTaxa)[counter])]] = PredictedPheno
      counter = counter +1
    }
    #How well do we predict observed lines at unobserved or observed timepoints:
    ObservedLinePrediction = FPCA$RecoveredCurves[c(1,28,55,83,123,205,300),] %>% as.data.frame() %>% 
      mutate(time = c(0,14,28,42,63,105,154)) %>% pivot_longer(cols = !time, names_to = 'taxa',values_to = 'GI3estimated') %>%
      merge(Truth, by = c('taxa','time')) %>% group_by(time,TP) %>% summarise(correlation = cor(GI3, GI3estimated)) %>% 
     mutate(fold = ii, Method = 'TrainingSetCorrelations')
    
    #GPP using predicted FPC's
    c(0,13.91,27.81,42.23,62.84,105.07,154)
    PredictedValues = as.matrix(FPCA$phi.fun.df[c(1,28,55,83,123,205,300),2:(numPcs+1)]) %*% t(as.matrix(PredParams[,-1])) +
      FPCA$EstimatedAndEmpiricalMu[c(1,28,55,83,123,205,300)+numObs,2]
    colnames(PredictedValues) = PredParams$taxa
    GPP = PredictedValues %>% data.frame() %>% mutate(time = c(0,14,28,42,63,105,154)) %>% pivot_longer(cols = !time, names_to = 'taxa',values_to = 'GI3estimated') %>%
      merge(Truth, by = c('taxa','time')) %>% group_by(time,TP) %>% summarise(correlation = cor(GI3, GI3estimated)) %>% 
      mutate(fold = ii,Method = 'GPP')
    
    #GPT using predicted time point values from the FPCA inputted lines
    GPT_training = t(FPCA$RecoveredCurves[c(1,28,55,83,123,205,300),]) %>% data.frame() %>% rename(TP1 = X1, TP2 = X2, TP3 = X3, TP4 = X4, TP5 = X5, TP6 = X6, TP7 = X7) 
    counter = 1
    GPTresults = data.frame(taxa = unique(GIlines$taxa)[testingnums])
    for (col in GPT_training){
        trained.Model = mixed.solve(col, Z = marker_train, K = NULL, SE =FALSE)
        PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u) + as.numeric(trained.Model$beta)
        print(trained.Model$beta)
        GPTresults[[(colnames(GPT_training)[counter])]] = PredictedPheno
        counter = counter +1
    }
    GPT = GPTresults %>% pivot_longer(!taxa, names_to = 'TP', values_to = 'GIestimate') %>% merge(.,Truth, by = c('taxa','TP')) %>%
      group_by(time,TP) %>% summarize(correlation = cor(GI3, GIestimate)) %>% mutate(fold = ii,Method = 'GPT')
    Results = rbind(Results, GPT, GPP,ObservedLinePrediction)
  }
  
  return(Results %>% mutate(Model = 'GIFPCA'))
  
}

# Drop one TP GI estimates FPCA #####
GPGIDropTPofObsFPCA =data.frame()
for (iii in c('None','TP2','TP3','TP4','TP5','TP6')){
  temp = GermFPCAcomplexGP(GIlines = all1920GIBlues.FPCAInput %>%  filter(taxa %in% CULines$taxa & taxa %nin% c('P_5_3','P_8_6'))%>%
                           filter(TP != iii) %>%
                           select(taxa,time,GI3),
                         myGD = myGD20_pruneCUfilt,
                         numfolds = 50,
                         datasplit = .80,
                         Truth = all1920GIBlues %>% select(taxa, time,TP, GI3)) 
  temp = temp %>% mutate(Dropped = iii)
  GPGIDropTPofObsFPCA = rbind(temp,GPGIDropTPofObsFPCA)
  
}

# Drop one TP GE estimates FPCA #####
GPGEDropTPofObsFPCA =data.frame()
for (iii in c('None','TP2','TP3','TP4','TP5','TP6')){
  temp =  GermFPCAcomplexGP(GIlines = All.BluesGE31920.FPCAInput %>%  filter(taxa %in% CULines$taxa & taxa %nin% c('P_5_3','P_8_6'))%>%
                            filter(TP != iii) %>%
                            select(taxa,time,GE3) %>% rename(GI3 = GE3),
                          myGD = myGD20_pruneCUfilt,
                          numfolds = 50,
                          datasplit = .80,
                          Truth = All.BluesGE31920 %>% select(taxa, time,TP, GE3) %>% rename(GI3 = GE3)) %>% mutate(Model = 'GEFPCA')
  
  temp = temp %>% mutate(Dropped = iii)
  GPGEDropTPofObsFPCA = rbind(temp,GPGEDropTPofObsFPCA)
}

# Save FPCA GP prediction outputs and plot them ####
FPCA_GEGI_predictions = rbind(GPGIDropTPofObsFPCA,
                              GPGEDropTPofObsFPCA) %>%
  mutate(ModelTrait = Model, Trait = substr(ModelTrait,1,2), Model = 'FPCA')

# setwd(rprojroot::find_rstudio_root_file())
# save(FPCA_GEGI_predictions, file = 'TimeSeriesAnalysis/Output/FPCA_GEGI_predictions.RData')

# Rectangle data frame for highlighting things #####
rectangles = FPCA_GEGI_predictions %>% select(TP, Trait, Dropped) %>% unique()%>%
  mutate(Top = 5, Bottom = -5, TraitLabel = paste0('Trait: ', Trait), TP = Dropped) %>%
  filter(Dropped !='None') %>% mutate(alpha = 0.05)

# GE Logistic drop a TP methods #######
# GE3 TP2 lines that fail
c("G_16_6","G_35_8","NY18101B_1",'NY18101B_2',"NY18102B_2",'NY18108B_3') # and others...
# GE3 TP3 Lines that fail
c("G_31_5","NY18102B_5","P_29_5","SB185_2")
# GE3 TP4 lines that fail
c("G_31_5",'P_16_1')
# GE3 TP5 lines that fail
c('P_29_5','P_32_5')
# GE3 TP6 lines that fail 
c("G_31_5")
GE3DropTPFittingFails = list(None = c("G_31_5","P_16_1",'P_29_5'),
  TP3 = c("G_31_5","NY18102B_5","P_29_5","SB185_2"),
                             TP4 = c("G_31_5",'P_16_1'),
                             TP5 = c('P_29_5','P_32_5'),
                             TP6 = c('G_31_5'))
# GE log complex drop one TP function #####
GELogComplexGPDropTP = function(GEinput = All.BluesGE31920 %>% filter(taxa %in% CULines$taxa &
                                                                        taxa %nin% c('P_5_3','P_8_6')) %>%
                                  filter(TP != iii & taxa %nin% GE3DropTPFittingFails[[iii]]),
                          myGD = myGD20_pruneCUfilt,
                          numfolds = 50,
                          datasplit = .80,
                          Truth = All.BluesGE31920 %>% select(taxa, time,TP, GE3)){
  
  GEParamsDropTP = GEinput %>% group_by(taxa) %>% 
    group_modify(~ broom::tidy(drm(GE3~time, data = .x,
                                   fct=LL.4(fixed = c(NA, NA, 1, NA),
                                            names = c('Rate','Lower','Upper','Centering'))))) %>%
    ungroup()
  
  
  GE3ParametersAndTPToCompare = GEParamsDropTP %>% select(taxa, term, estimate) %>% pivot_wider(values_from = 'estimate',names_from = 'term') %>%
    select(taxa, Rate, Lower, Centering) %>%
    mutate(TP1 = Lower+(1-Lower)/(1+exp(Rate*(log(0)-log(Centering)))),
           TP2 = Lower+(1-Lower)/(1+exp(Rate*(log(14)-log(Centering)))),
           TP3 = Lower+(1-Lower)/(1+exp(Rate*(log(28)-log(Centering)))),
           TP4 = Lower+(1-Lower)/(1+exp(Rate*(log(42)-log(Centering)))),
           TP5 = Lower+(1-Lower)/(1+exp(Rate*(log(63)-log(Centering)))),
           TP6 = Lower+(1-Lower)/(1+exp(Rate*(log(105)-log(Centering)))),
           TP7 = Lower+(1-Lower)/(1+exp(Rate*(log(154)-log(Centering)))))
  
  GEinput = GEinput %>% arrange(taxa, time)
  checkorder = sum(GEinput %>% select(taxa) %>% unique() == 
                     as.character(((myGD %>% filter(taxa %in% GEinput$taxa) %>% arrange(taxa))$taxa)))
  myGD = (myGD %>% filter(taxa %in% GEinput$taxa) %>% arrange(taxa))[,-1]-1
  Results = data.frame()
  set.seed(1)
  num_entries = dim(GE3ParametersAndTPToCompare)[1]
  Training_size = round(datasplit*num_entries,0)

  if(checkorder != num_entries){
    stop('myGD and the phenotypes are not aligning correctly, ensure order')
  }
  
  # Lets start with things
   for (ii in 1:numfolds) {
    start_time <- Sys.time()
    trainingnums = sort(sample(1:num_entries,Training_size))
    testingnums = setdiff(1:num_entries,trainingnums)
    
    trainingSet = GE3ParametersAndTPToCompare[trainingnums,]
    testingSet = GE3ParametersAndTPToCompare[testingnums,]
    marker_train = myGD[trainingnums,]
    marker_test = myGD[testingnums,]
    PredParams = data.frame(taxa = testingSet$taxa)
    counter = 2
    for (col in trainingSet[,-1]){
      trained.Model = mixed.solve(col, Z = marker_train, K = NULL, SE =FALSE)
      PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u) + as.numeric(trained.Model$beta)
      PredParams[[(colnames(trainingSet)[counter])]] = PredictedPheno
      counter = counter +1
    }
    #LEts look here at the per TP estimates using the timepoint estimates
    GPT = PredParams %>% select(!c(Rate, Lower, Centering))%>% 
      pivot_longer(cols = !taxa, names_to = 'TP', values_to = 'GE3perTPestimate')%>%
      merge(.,Truth, by = c('taxa','TP')) %>% group_by(TP,time) %>% 
      summarise(correlation = cor(GE3perTPestimate,GE3, use = 'complete.obs')) %>% mutate(fold = ii,Method = 'GPT')
    
    #Here we calculate the correlation between the parameter estimated GE using GP'ed parameter estimates
    PredParams.l = PredParams %>% select(taxa, Rate, Lower, Centering)%>% 
      pivot_longer(cols = !taxa, names_to = 'term', values_to = 'estimate')
    # 0  14  28  42  63 105 154 are the times we need
    GEperTPestimates = data.frame()
    for (i in unique(PredParams.l$taxa)){
      time = c(0,14,  28,  42,  63, 105, 154)
      tmp = PredParams.l %>% filter(taxa ==i)
      y = tmp$estimate[2]+(1-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[3]))))
      GEperTPestimates = rbind(GEperTPestimates, data.frame(taxa = i, time = time, GE3estimate = y))
    }
    GPP = GEperTPestimates %>% merge(.,
                                     Truth, 
                                     by = c('taxa','time')) %>% group_by(TP,time) %>% 
      summarise(correlation = cor(GE3estimate,GE3, use = 'complete.obs'))%>%mutate(fold = ii,Method = 'GPP')
    ObsInObservedTimeCor = trainingSet %>% select(taxa,TP1,TP2,TP3,TP4,TP5,TP6,TP7) %>% 
      pivot_longer(cols = starts_with('TP'),names_to = 'TP',values_to = 'GEestimate') %>%
      merge(.,
            Truth, 
            by = c('taxa','TP')) %>% group_by(TP,time) %>% 
      summarise(correlation = cor(GEestimate,GE3))%>%mutate(fold = ii,Method = 'TrainingSetCorrelations')
    
    Results = Results %>% rbind(.,GPT,GPP,ObsInObservedTimeCor)
    end_time <- Sys.time()
    print(end_time-start_time)
  }
  return(Results %>% mutate(Model = 'GELogistic'))
}

# GI logistic drop a TP and full models ############
# Taxa that make the models fail 
# TP3
c("G_31_9","NY18113B_2","NY18115B_2","NY18124B_2","NY18124B_4",
  "P_14_6","P_2_5","P_26_3","P_28_5","P_32_1","P_36_6","SB184R_4",
  "SP333R_1","SP353R_2","SR6125_3","ST1444R_1")
#TP4
c("P_26_2","P_36_6","SB182R_3","SN341_3")
# TP5
c("G_31_5","NY18120B_2","NY18120B_3","NY18120B_4","NY18125B_3",
  "P_2_5", "P_26_3","SB193R_1","SB581R_3","SN813_4","SP575_1","ST1431R_1")
# TP6
c("G_31_5","NY18120B_2","P_26_3","SB182R_3","SG5123_1","SG592R_1","SP575_1")
GI3DropTPFittingFails = list(None = c('G_31_5', 'Megs_song','NY18120B_4','NY18125B_1','Oderbrucker',
                             'P_2_5','P_26_3','P_36_6', 'SB193R_1','SG514R_4','SN873_3',
                             'ST1431R_1'),
  TP3 = c("G_31_9","NY18113B_2","NY18115B_2","NY18124B_2","NY18124B_4",
          "P_14_6","P_2_5","P_26_3","P_28_5","P_32_1","P_36_6","SB184R_4",
          "SP333R_1","SP353R_2","SR6125_3","ST1444R_1"),
  TP4 = c("P_26_2","P_36_6","SB182R_3","SN341_3"),
  TP5 = c("G_31_5","NY18120B_2","NY18120B_3","NY18120B_4","NY18125B_3",
          "P_2_5", "P_26_3","SB193R_1","SB581R_3","SN813_4","SP575_1","ST1431R_1"),
  TP6 = c("G_31_5","NY18120B_2","P_26_3","SB182R_3","SG5123_1","SG592R_1","SP575_1"))

# GI logistic drop one TP function for calling######
GILogComplexGPdropTP = function(GIinput = all1920GIBlues %>%select(taxa, time,TP, GI3) %>%
                                  filter(taxa %in% CULines$taxa & taxa %nin% c('P_5_3','P_8_6')) %>%
                                  filter(TP != iii & taxa %nin% GI3DropTPFittingFails[[iii]]),
                          myGD = myGD20_pruneCUfilt,
                          numfolds = 50,
                          datasplit = .80,
                          Truth = all1920GIBlues %>% select(taxa, time,TP, GI3)){

  
  GIParamsDropTP = GIinput %>% group_by(taxa) %>% 
    group_modify(~ broom::tidy(drm(GI3~time, data = .x,
                                   fct=LL.4(fixed = c(NA, NA, NA, NA),
                                            names = c('Rate','Lower','Upper','Centering'))))) %>%
    ungroup()
  

  GI3ParametersAndTPToCompare = GIParamsDropTP %>% select(taxa, term, estimate) %>%
    pivot_wider(names_from = 'term',values_from = 'estimate') %>%
                            mutate(TP1 = Lower+(Upper-Lower)/(1+exp(Rate*(log(0)-log(Centering)))),
                                   TP2 = Lower+(Upper-Lower)/(1+exp(Rate*(log(14)-log(Centering)))),
                                   TP3 = Lower+(Upper-Lower)/(1+exp(Rate*(log(28)-log(Centering)))),
                                   TP4 = Lower+(Upper-Lower)/(1+exp(Rate*(log(42)-log(Centering)))),
                                   TP5 = Lower+(Upper-Lower)/(1+exp(Rate*(log(63)-log(Centering)))),
                                   TP6 = Lower+(Upper-Lower)/(1+exp(Rate*(log(105)-log(Centering)))),
                                   TP7 = Lower+(Upper-Lower)/(1+exp(Rate*(log(154)-log(Centering)))))# Lets start with things
  Results = data.frame()
  set.seed(1)
  num_entries = dim(GI3ParametersAndTPToCompare)[1]
  Training_size = round(datasplit*num_entries,0)

  GIinput = GIinput %>% arrange(taxa, time)
  checkorder = sum(GIinput %>% select(taxa) %>% unique() == 
                     as.character(((myGD %>% filter(taxa %in% GIinput$taxa) %>% arrange(taxa))$taxa)))
  myGD = (myGD %>% filter(taxa %in% GIinput$taxa) %>% arrange(taxa))[,-1]-1
  
  if(checkorder != num_entries){
    stop('myGD and the phenotypes are not aligning correctly, ensure order')
  }
  
  for (ii in 1:numfolds) {
    start_time <- Sys.time()
    trainingnums = sort(sample(1:num_entries,Training_size))
    testingnums = setdiff(1:num_entries,trainingnums)

    trainingSet = GI3ParametersAndTPToCompare[trainingnums,]
    testingSet = GI3ParametersAndTPToCompare[testingnums,]
    marker_train = myGD[trainingnums,]
    marker_test = myGD[testingnums,]
    PredParams = data.frame(taxa = testingSet$taxa)
    counter = 2
    for (col in trainingSet[,-1]){
      trained.Model = mixed.solve(col, Z = marker_train, K = NULL, SE =FALSE)
      PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u) + as.numeric(trained.Model$beta)
      PredParams[[(colnames(trainingSet)[counter])]] = PredictedPheno
      counter = counter +1
    }
    #LEts look here at the per TP estimates using the timepoint estimates
    GPT = PredParams %>% select(!c(Rate, Lower, Centering, Upper))%>%
      pivot_longer(cols = !taxa, names_to = 'TP', values_to = 'GI3perTPestimate')%>%
      merge(.,Truth, by = c('taxa','TP')) %>% group_by(TP,time) %>%
      summarise(correlation = cor(GI3perTPestimate,GI3, use = 'complete.obs')) %>% mutate(fold = ii,Method = 'GPT')

    #Here we calculate the correlation between the parameter estimated GE using GP'ed parameter estimates
    PredParams.l = PredParams %>% select(taxa, Rate, Lower, Upper, Centering)%>%
      pivot_longer(cols = !taxa, names_to = 'term', values_to = 'estimate')
    # 0  14  28  42  63 105 154 are the times we need
    GIperTPestimates = data.frame()
    for (i in unique(PredParams.l$taxa)){
      time = c(0,14,  28,  42,  63, 105, 154)
      tmp = PredParams.l %>% filter(taxa ==i)
      y = tmp$estimate[2]+(tmp$estimate[3]-tmp$estimate[2])/(1+exp(tmp$estimate[1]*(log(time)-log(tmp$estimate[4]))))
      GIperTPestimates = rbind(GIperTPestimates, data.frame(taxa = i, time = time, GI3estimate = y))
    }
    GPP = GIperTPestimates %>% merge(.,
                                     Truth,
                                     by = c('taxa','time')) %>% group_by(TP,time) %>%
      summarise(correlation = cor(GI3estimate,GI3, use = 'complete.obs'))%>%mutate(fold = ii,Method = 'GPP')
    #observed lines accuracys
    ObsInObservedTimeCor = trainingSet %>% select(taxa,TP1,TP2,TP3,TP4,TP5,TP6,TP7) %>% 
      pivot_longer(cols = starts_with('TP'),names_to = 'TP',values_to = 'GIestimate') %>%
      merge(.,
            Truth, 
            by = c('taxa','TP')) %>% group_by(TP,time) %>% 
      summarise(correlation = cor(GIestimate,GI3, use = 'complete.obs'))%>%
      mutate(fold = ii,Method = 'TrainingSetCorrelations')
    
    
    Results = Results %>% rbind(.,GPT,GPP,ObsInObservedTimeCor)
    end_time <- Sys.time()
    print(end_time-start_time)
  }
  return(Results %>% mutate(Model = 'GILogistic'))
}
# Use the function to calculate the Drop one TP results for GI logistic. #####
GPGIDropTPofObsLog = data.frame()
for (iii in c('None','TP3','TP4','TP5','TP6')){
  temp = GILogComplexGPdropTP(GIinput = all1920GIBlues %>%select(taxa, time,TP, GI3) %>%
                                filter(taxa %in% CULines$taxa & taxa %nin% c('P_5_3','P_8_6')) %>%
                                filter(TP != iii & taxa %nin% GI3DropTPFittingFails[[iii]]),
                              myGD = myGD20_pruneCUfilt,
                              numfolds = 50,
                              datasplit = .80,
                              Truth = all1920GIBlues %>% select(taxa, time,TP, GI3))
  temp = temp %>% mutate(Dropped = iii)
  GPGIDropTPofObsLog = rbind(temp,GPGIDropTPofObsLog)
}

# Use the function to calculate the Drop one TP results for GE logistic. #####
GPGEDropTPofObsLog = data.frame()
for (iii in c('None','TP3','TP4','TP5','TP6')){
  temp = GELogComplexGPDropTP(GEinput = All.BluesGE31920 %>% filter(taxa %in% CULines$taxa &
                                                                      taxa %nin% c('P_5_3','P_8_6')) %>%
                                filter(TP != iii & taxa %nin% GE3DropTPFittingFails[[iii]]),
                              myGD = myGD20_pruneCUfilt,
                              numfolds = 50,
                              datasplit = .80,
                              Truth = All.BluesGE31920 %>% select(taxa, time,TP, GE3))
  temp = temp %>% mutate(Dropped = iii)
  GPGEDropTPofObsLog = rbind(temp,GPGEDropTPofObsLog)
}

# Bind the Ge and GI logistic together and save them! #####
Logistic_GEGI_predictions = rbind(GPGEDropTPofObsLog,
                                  GPGIDropTPofObsLog) %>%
  mutate(ModelTrait = Model, Trait = substr(ModelTrait,1,2), Model = 'Logistic')

# setwd(rprojroot::find_rstudio_root_file())
# save(Logistic_GEGI_predictions, file = 'TimeSeriesAnalysis/Output/Logistic_GEGI_predictions.RData')

# Two separate grpahs for logistic and FPCa then patchwork together #####
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
# setwd('SpringBarley/Analysis/All_taxa/FinalOutputpngs/')

# png('FPCA_LOG_GEGI_GP.png', 800,1000, res = 120)
FPCA_GEGI_predictions %>%
  mutate(MethodLabel = mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                 to = c('Training Set\nAccuracy',
                                        'Testing Set\nGPP Accuracy', 'Testing Set\nGPT Accuracy')),
         TraitLabel = paste0('Model: FPCA; Trait: ', Trait))%>%
  ggplot(aes(x = TP, y = correlation, color = MethodLabel, fill = MethodLabel))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(TraitLabel),rows = vars(Dropped)) +
  geom_hline(yintercept = 0, color = 'black') +
  ylab('Accuracy') +xlab('TP of Prediction')+
  labs(color = 'Prediction \nMethod')+ guides(fill = FALSE) +
  geom_tile(data = rectangles %>% mutate(TraitLabel = paste0('Model: FPCA; Trait: ', Trait)),
            inherit.aes = FALSE,
            aes(x = TP, fill = NA,y = 0), alpha = 0.02,
            height = 2, width = 1, color = 'lightgrey') +
  theme_bw(base_size = 9)+
Logistic_GEGI_predictions%>%
  mutate(MethodLabel = mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                 to = c('Training Set\nAccuracy',
                                        'Testing Set\nGPP Accuracy', 'Testing Set\nGPT Accuracy')),
         TraitLabel = paste0('Model: Logistic; Trait: ', Trait))%>%
  ggplot(aes(x = TP, y = correlation, color = MethodLabel, fill = MethodLabel))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(TraitLabel),rows = vars(Dropped)) +
  geom_hline(yintercept = 0, color = 'black') +
  ylab('Accuracy') +xlab('TP of Prediction')+
  labs(color = 'Prediction \nMethod')+ guides(fill = FALSE) +
  geom_tile(data = rectangles %>% filter(Dropped != 'TP2')%>% 
              mutate(TraitLabel = paste0('Model: Logistic; Trait: ', Trait)),
            inherit.aes = FALSE,
            aes(x = TP, fill = NA,y = 0), alpha = 0.02,
            height = 2, width = 1, color = 'lightgrey') +
  theme_bw(base_size = 9) +
  plot_layout(design = 'A\nB', guides = 'collect')
# dev.off()

# Model:Trait across the Top with Trait Ordering #####
# jpeg('FPCA_LOG_GEGI_GPv2.jpg', 1000,800, res = 120)
FPCALogisticPrediction %>%
  ggplot(aes(x = TP, y = correlation, color = MethodLabel))+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a'))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(TraitLabel),rows = vars(Dropped)) +
  geom_hline(yintercept = 0, color = 'black') +
  ylab('Accuracy') +xlab('TP of Prediction')+
  labs(color = 'Prediction \nMethod') +
  geom_tile(data = Rectangles2,
            inherit.aes = FALSE,
            aes(x = TP, fill = NA,y = 0), alpha = 0.05,
            height = 2, width = 1, fill = 'lightgrey') +
  theme_bw(base_size = 9) +
  ylab('Prediction accuracy')
# dev.off()


library(ggh4x)
# jpeg('FPCA_LOG_GEGI_GPv3facet_strip.jpg', 1000,800, res = 120)
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
# dev.off()

FPCALogisticPrediction %>% mutate('TP Dropped' = 'Time point masked') %>%
  ggplot(aes(x = TP, y = correlation, color = MethodLabel))+
  geom_boxplot(outlier.size = 0.5) +
  facet_nested(`TP Dropped`+Dropped  ~ TraitLabel, scales = 'free_y')+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a'))+
  geom_hline(yintercept = 0, color = 'black') +
  xlab('Time point of prediction')+
  labs(color = 'Prediction \nMethod') +
  # geom_tile(data = Rectangles2 %>% mutate('TP Dropped' = 'Time point masked'),
  #           inherit.aes = FALSE,
  #           aes(x = TP, fill = NA,y = 0), alpha = 0.05,
  #           height = 2, width = 1, fill = 'lightgrey') +
  theme_bw(base_size = 10) +
  ylab('Prediction accuracy')





# plot by TP:Model x- axis #####
# png('FPCA_LOG_GEGI_GPv3.png', 1000,800, res = 120)
FPCALogisticPrediction%>%
  ggplot(aes(x = TPModel, y = correlation, color = MethodLabel))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(Trait),rows = vars(Dropped)) +
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a'))+
  geom_hline(yintercept = 0, color = 'black') +
  ylab('Accuracy') +xlab('TP of Prediction/Model')+
  labs(color = 'Prediction \nMethod') +
  geom_tile(data = Rectangles2,
            inherit.aes = FALSE,
            aes(x = TPModel, fill = NA,y = 0), alpha = 0.05,
            height = 2, width = 1, fill = 'lightgrey') +
  theme_bw(base_size = 9) +theme(axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5))
# dev.off()
######

Logistic_GEGI_predictions %>%
  mutate(TraitLabel = paste0('Model: Logistic; Trait: ', Trait)) %>%
  rbind(FPCA_GEGI_predictions %>%
          mutate(TraitLabel = paste0('Model: FPCA; Trait: ', Trait))) %>% 
  mutate(TraitLabel = factor(TraitLabel,levels = c('Model: FPCA; Trait: GE',
                                                   'Model: Logistic; Trait: GE',
                                                   'Model: FPCA; Trait: GI',
                                                   'Model: Logistic; Trait: GI')),
         MethodLabel = mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                 to = c('Training Set\nAccuracy',
                                        'Testing Set\nGPP Accuracy', 'Testing Set\nGPT Accuracy')),
         TPModel = paste0(TP,':',Model)) %>%
  group_by(Model, Trait, Dropped, TP, Method) %>%
  summarise(PA = mean(correlation), stand.dev = sd(correlation), Min = min(correlation), Max = max(correlation)) %>%
  mutate(Line_Status = ifelse(Method %in% c('GPT','GPP'), 'Unobs','Obs'),
         Timepoint_status = ifelse(TP==Dropped, 'Unobs','Obs'),
         Line_time_status = paste0(Line_Status,'-',Timepoint_status)) 

a =Logistic_GEGI_predictions%>%
  mutate(MethodLabel = mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                 to = c('Training Set PA',
                                        'Testing Set GPP PA', 'Testing Set GPT PA')),
         TraitLabel = paste0('Model: Logistic; Trait: ', Trait)) %>%
  rbind(FPCA_GEGI_predictions %>%
          mutate(MethodLabel = mapvalues(Method, from = c('TrainingSetCorrelations','GPP','GPT'), 
                                         to = c('Training Set PA',
                                                'Testing Set GPP PA', 'Testing Set GPT PA')),
                 TraitLabel = paste0('Model: FPCA; Trait: ', Trait))) %>%
  group_by(Model, Trait, Dropped, TP, MethodLabel) %>%
  summarise(PA = round(mean(correlation),2)) %>% select(TP,Trait, Model,MethodLabel, Dropped, PA) %>%
  pivot_wider(id_cols = c('TP','Dropped'),
              names_from = c('Trait','Model','MethodLabel'), names_sep = ' ', values_from = 'PA') %>%
  rename('TP masked'=Dropped) 

# write.table(a, file = 'SpringBarley/Manuscripts/PA_Summary.csv',sep = ',')
FPCALogisticPrediction %>%
  mutate(Line_Status = ifelse(Method %in% c('GPT','GPP'), 'Unobs','Obs'),
         Timepoint_status = ifelse(TP==Dropped, 'Unobs','Obs'),
         Line_time_status = paste0(Line_Status,'-',Timepoint_status)) %>%
  filter(Method != 'GPT') %>%
  ggplot(aes(x = TP, y = correlation, color = Line_time_status, shape = Model))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(TraitLabel),rows = vars(Dropped)) +
  geom_hline(yintercept = 0, color = 'black') +
  ylab('Accuracy') +xlab('TP of Prediction')+
  labs(color = 'Prediction \nMethod')
