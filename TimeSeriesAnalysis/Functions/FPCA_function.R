# This function perfrom FPCA on curves given some inputs:

FPCA_function = function(dfTaxaTraitTime,  # dfTaxaTraitTime is a df with taxa, trait, and time. Columns must be named 'taxa' and 'time'.
                         # !!! FILTER YOUR DATA BEFORE ENTERING IT INTO THIS FUNCTION, DATA MUST BE BALANCED !!!
                         Trait, #Trait name must be entered as a character ie Trait = 'GE3'
                         NumKnots, # NumKnots is the number of interior knots to be fitted
                         order, # Order is the degree of the polynomial to be fit to the data. 
                         NumObsevationsPerLine = dim(dfTaxaTraitTime %>% filter(taxa=='AAC_Synergy'))[1]){
  
  dfTaxaTraitTime  = dfTaxaTraitTime %>% arrange(taxa, time)
  NumLines = dim(dfTaxaTraitTime %>% group_by(taxa) %>% dplyr::summarize())[1] 
  
  trait.Y.vec = dfTaxaTraitTime[[Trait]]	# vectorize the observation values
  time.vec = dfTaxaTraitTime$time	# vectorize the time points
  
  ylimit = range(trait.Y.vec)	# range of observation values
  tlimit = range(time.vec)	# range of measurement time points
  ## draw a picture for the raw data and empirical mean curve #######
  time.points = unique(sort(time.vec))
  mu.emp = rep(0,length(time.points))
  for (i in 1:length(time.points)){
    mu.emp[i] = mean(trait.Y.vec[time.vec==time.points[i]])
  }
  print(time.points)
  print(mu.emp)
  # raw data with the empirical mean function
 EmpiricalMeanplot =  ggplot(dfTaxaTraitTime)+geom_line(aes(x = time, y = .data[[Trait]], group = taxa))+
    geom_line(data = data.frame(time = time.points, muEmp = mu.emp), aes(x = time, y = muEmp), color = 'green', size = 2)
  #### Find estimated mean function ####
  Knot.int = NumKnots	# Knot.int is the number of interior knots
  knots.loc = tlimit[1] + (tlimit[2]-tlimit[1])*(1:Knot.int)/(1+Knot.int)	# time points of interios knots
  K_num_basis = length(knots.loc) + order 				# number of basis functions
  basis.fnct = create.bspline.basis(tlimit,K_num_basis,norder=order)

  print(paste('Interior knot at ', knots.loc))
  
  Omega = inprod(basis.fnct,basis.fnct,2,2)			# penalty matrix
  inte1 = kronecker(inprod(basis.fnct,basis.fnct,2,2),inprod(basis.fnct,basis.fnct))
  inte2 = kronecker(inprod(basis.fnct,basis.fnct,1,1),inprod(basis.fnct,basis.fnct,1,1))
  inte3 = kronecker(inprod(basis.fnct,basis.fnct),inprod(basis.fnct,basis.fnct,2,2))
  Omega2 = inte1+2*inte2+inte3
  
  N.obs = length(trait.Y.vec)
  Xmat = matrix(0,N.obs,K_num_basis)
  start.temp = 1
  for (i in 1:NumLines){
    n.i = length(time.points)
    Xmat[(start.temp:(start.temp+n.i-1)),] = eval.basis(time.points,basis.fnct,0)
    start.temp = start.temp+n.i
  }
  
  ### Penalized least squares estimates ###
  lam = tuning_nointer(-10,15,Omega,Xmat,trait.Y.vec)	# tunning parameter
  # print(log(lam))
  bhat = solve(t(Xmat)%*%Xmat+adiag(Omega*lam))%*%t(Xmat)%*%trait.Y.vec
  print('Estimating the mean function...')
  J = 300
  tt = seq(tlimit[1], tlimit[2], length.out=J)	# evaluation time points
  BS = eval.basis(tt,basis.fnct,0)		# evaluation of the basis functions over the evaluation time points
  mu = BS%*%bhat	# estimate mean function
  
  EstimatedAndEmpiricalMu = data.frame(time = as.vector(x = c(time.points,tt)),
                                Mu = as.vector(c(mu.emp, mu)),
                                EstimationMethod = c(rep('Empirical',length(time.points)),rep('Fitted',length(tt))))
  EstimatedMeanPlot = ggplot(dfTaxaTraitTime)+geom_line(aes(x = time, y = .data[[Trait]], group = taxa))+
    geom_line(data = EstimatedAndEmpiricalMu, 
              aes(x = time, y = Mu, group = EstimationMethod, color = EstimationMethod), size = 2)
  #### Estimate the PCA  eigenfucntions from the mean function ####
  print('estimating eigenFunctions...')
  tempPCA_funtionsOutput = pca_fun(trait.Y.vec, T.vec = time.vec,
                                   Xmat, bhat, K_num_basis, NumLines,
                                   BS, basis.fnct, Omega, Omega2,
                                   NumObsPerLine = NumObsevationsPerLine,J)
  
  v1 = tempPCA_funtionsOutput[[1]]
  V1 = tempPCA_funtionsOutput[[2]]
  sigma.e2.hat = tempPCA_funtionsOutput[[3]]
  print('Esitmated... Moving to FPCs per line')
  
  # The eigen-function matrix
  Phi1 = list(NumLines)
  start_temp = 1
  for (i in unique(dfTaxaTraitTime$taxa)){
    temp.taxa = dfTaxaTraitTime %>% filter(taxa == i) %>% dplyr::select(time)
    Phi1[[start_temp]] =  eval.basis(temp.taxa$time,basis.fnct,0)%*%V1
    start_temp= start_temp+1
  }
  
  # Choose the number of PCs by PVE method
  v1.trim = v1[v1>0]
  n1 = length(v1.trim)
  s1 = sum(v1.trim)
  
  phi.fun = BS%*%V1
  phi.fun.df = data.frame(time = tt, phi.fun) 
  colnames(phi.fun.df)[2:(dim(phi.fun)[2]+1)] = gsub("X", "Phi", colnames(phi.fun.df)[2:(dim(phi.fun)[2]+1)])
  
  phi.fun.plot = phi.fun.df %>% pivot_longer(cols = starts_with("Phi")) %>% 
    ggplot(aes(x = time, y = value, color = name, group = name))+geom_line()+
    labs(color = 'Phi \n Function')
  phi.funwithMu.plot = phi.fun.df %>% pivot_longer(cols = starts_with('Phi')) %>%  arrange(name) %>%
    mutate(value = value + rep(mu,dim(phi.fun)[2])) %>% 
    ggplot(aes(x = time, y = value, color = name, group = name))+geom_line() + labs(color = 'Phi \n Function')
  
  print('Phi functions plotted')
  print('Predicting PCs for each line... this may take a minute')
  
  PCs = pca_score(N = N.obs,Y.vec = trait.Y.vec, Xmat,bhat,
                       K.est = dim(phi.fun)[2],v1,Phi1,NumLines,
                       sigma.e2.hat,
                       NumObsPerLine = NumObsevationsPerLine)
  summary(PCs)
  
  
  PCs_withTaxa = data.frame(taxa = unique(dfTaxaTraitTime$taxa), PCs)
  colnums = dim(phi.fun)[2]
  colnames(PCs_withTaxa)[2:(colnums+1)] = gsub("X", "FPC", colnames(PCs_withTaxa)[2:(colnums+1)])
  head(PCs_withTaxa)
  
  colindex = length(v1.trim)+1

  MuFunction = EstimatedAndEmpiricalMu %>% dplyr::filter(EstimationMethod == 'Fitted')
  # Lets just take the first 3 FPCs to reconstruct the GI3 Curves. 
  RecoveredCurves  = as.matrix(phi.fun.df[,2:colindex]) %*% t(as.matrix(PCs_withTaxa[,2:colindex])) + 
    matrix(rep(MuFunction$Mu,dim(PCs_withTaxa)[1]),ncol = dim(PCs_withTaxa)[1])
  colnames(RecoveredCurves) = PCs_withTaxa$taxa
  RecoveredCurvePlot = as.data.frame(RecoveredCurves) %>% mutate(time = phi.fun.df$time) %>% 
    pivot_longer(cols = 1:dim(PCs_withTaxa)[1], names_to = 'taxa') %>%
    ggplot()+geom_line(aes( x= time,y = value, group =taxa))+
    geom_point(data = dfTaxaTraitTime, aes(x=time,y= .data[[Trait]], group = taxa),color = 'blue')
  
  
  
  
  
  return(list('EstimatedAndEmpiricalMu' = EstimatedAndEmpiricalMu,
              'phi.fun.df' = phi.fun.df, 
              'PCs_withTaxa' =  PCs_withTaxa, 
              'v1' = tempPCA_funtionsOutput[[1]],
              'V1' = tempPCA_funtionsOutput[[2]],
              'sigma.e2.hat' = tempPCA_funtionsOutput[[3]],
              'EmpiricalMeanplot' = EmpiricalMeanplot,
              'EstimatedMeanPlot'= EstimatedMeanPlot,
              'phi.fun.plot' = phi.fun.plot,
              'phi.funwithMu.plot' = phi.funwithMu.plot,
              'RecoveredCurves' = RecoveredCurves,
              'RecoveredCurvePlot' = RecoveredCurvePlot))
  
}



