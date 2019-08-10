
# y = numeric vector, response variable
# gen = genotypic matrix, preferably coded as 012
# k = integer, number of folds for cross-validation
# Seeds = numeric vector, random seed to each cross-validation
# SET = factor/vector, if provided CVs are performed as leave-level-out
# box = if FALSE, the function only returns the average across CVs

CV = function(y,gen,k=5,Seeds=1:20,SET=NULL,box=TRUE){

  # Remove missing
  if(anyNA(y)){
    wna = which(is.na(y))
    y = y[-wna]
    gen = gen[-wna,]
    if(!is.null(SET)) SET = SET[-wna]
  }
  
  # Kernels
  timeE2 = system.time(E2 <- as.matrix(dist(gen)^2))[[3]]
  timeK = system.time(K <- exp(-(E2/mean(E2))))[[3]]
  timeEK = system.time(eK <- eigen(K,symmetric=T))[[3]]
  timeDecomp = timeE2+timeK+timeEK
  
  timeCNT = system.time(Z <- apply(gen,2,function(x) x-mean(x)))[[3]]
  timeZZ = system.time(G <- tcrossprod(Z))[[3]]
  timeSca = system.time(G <- G/(2*sum(apply(Z,2,var))))[[3]]
  timeBuildG = timeCNT+timeZZ+timeSca
  
  Y = y
  N = nrow(gen)
  cvs = length(Seeds)
  cat('Done with initial settings\n')
  
  # Cross-validation function
  folds = function(Seeds,SET=NULL){
    
    # Loading packages
    require(BGLR,quietly = T)
    require(ranger,quietly = T)
    require(glmnet,quietly = T)
    require(gbm,quietly = T)
    require(rrBLUP,quietly = T)
    require(bWGR,quietly = T)
    require(pls,quietly = T)
    require(kernlab,quietly = T)
    require(VIGoR,quietly = T)
    require(EBglmnet,quietly = T)
    
    Rcpp::cppFunction('
                      SEXP FLM(NumericVector y, NumericMatrix X){
                      // Convergence settings
                      int maxit = 300; double tol = 10e-8;
                      // Initial settings and starting values
                      int p=X.ncol(), n=X.nrow(), numit=0;
                      double b0,b1,eM,Ve,cnv=1,mu=mean(y),Lmb2=0;
                      NumericVector e=y-mu,Vb(p),b(p),fit(n);
                      // Cross-products and shape parameter
                      NumericVector xx(p),sx(p),bc(p);
                      for(int k=0; k<p; k++){
                      xx[k]=sum(X(_,k)*X(_,k));
                      if(xx[k]==0) xx[k]=0.1;
                      Lmb2=Lmb2+var(X(_,k));}
                      NumericVector iTau2=p+Lmb2;
                      // Looping across parameters until convergence
                      while(numit<maxit){   
                      // Updating markers effects
                      bc=b+0; for(int j=0; j<p; j++){ b0=b[j];
                      b1=(sum(X(_,j)*e)+xx[j]*b0)/(iTau2(j)+xx(j));
                      b[j]=b1; e=e-X(_,j)*(b1-b0);}
                      // Updating intercept
                      eM=mean(e); mu=mu+eM; e=e-eM;
                      // Updating variance components
                      Ve=sum(e*y)/(n-1); 
                      Vb=b*b+Ve/(xx+iTau2);
                      iTau2=sqrt(Lmb2*Ve/Vb);
                      // Check parameters convergence
                      ++numit; cnv=sum(abs(bc-b));
                      if(cnv<tol){break;}}
                      // Fit model
                      for(int k=0; k<n; k++){fit[k]=mu+sum(X(k,_)*b);}
                      // Return output
                      return List::create(Named("mu")=mu,Named("b")=b);}')
    
    # Begin folds / Pre-defined subset
    if(is.null(SET)){
      set.seed(Seeds)
      Nk = round(N/k)
      w = sample(1:N,Nk)
      y[w] = NA
    }else{
      w = SET
      y[w] = NA
    }
    
    # RF
    timeQ_RndFrst = system.time(ff11_RndFrst <- predict(ranger(y~.,data.frame(y=y[-w],gen[-w,]),importance='impurity',num.trees=1000),data.frame(gen[w,]))$predictions)[[3]]
    cat('ranger\n')
    
    # Kernels
    timeQ_EpsSVR = system.time( f22b <- predict(ksvm(gen,y,type= "eps-svr"),gen)[,1])[[3]]; ff22b_EpsSVR = f22b[w]
    timeQ_NuSVR = system.time( f22c <- predict(ksvm(gen,y,type= "nu-svr"),gen)[,1])[[3]]; ff22c_NuSVR = f22c[w]
    cat('kernlab\n')
    
    # Bayesian
    timeQ_BayesA = system.time( ff3b_BayesA <- BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesA')),verbose=F)$yHat[w] )[[3]]
    timeQ_BayesB = system.time( ff2b_BayesB <- BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesB')),verbose=F)$yHat[w] )[[3]]
    timeQ_BayesC = system.time( ff2c_BayesC <- BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesC')),verbose=F)$yHat[w] )[[3]]
    timeQ_BL = system.time( ff4b_BL <- BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BL',type='beta')),verbose=F)$yHat[w])[[3]]
    timeQ_BRR = system.time( ff5b_BRR <- BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BRR')),verbose=F)$yHat[w])[[3]]
    timeQ_RKHS = system.time( ff3_RKHS <- BGLR(rmExistingFiles = F,y,ETA=list(list(V=eK$vectors,d=eK$values,model='RKHS')),verbose=F)$yHat[w] )[[3]] + timeDecomp
    file.remove(list.files(pattern='dat'))
    cat('BGLR\n')
    
    # Bayesian 2
    timeQ_BayesCpi = system.time( ff4_BayesCpi <- BayesCpi( y[-w], gen[-w,] ))[[3]]; ff4_BayesCpi = gen[w,] %*% ff4_BayesCpi$b 
    timeQ_BayesDpi = system.time( ff4_BayesDpi <- BayesDpi( y[-w], gen[-w,] ))[[3]]; ff4_BayesDpi = gen[w,] %*% ff4_BayesDpi$b
    cat('bWGR pi\n')
    
    # L1 L2
    timeCV1 = system.time( cv1 <- cv.glmnet(x=gen[-w,],y=y[-w],alpha=1) )[[3]]
    lmb1 <- cv1$lambda.min;
    timeFIT1 = system.time( f4 <- glmnet(x=gen[-w,],y=y[-w],lambda = lmb1,alpha = 1) )[[3]]
    ffw4_LASSO <- c(predict(f4,gen[w,]))
    timeQ_LASSO = timeCV1 + timeFIT1
    #
    timeCV2 = system.time( cv2 <- cv.glmnet(x=gen[-w,],y=y[-w],alpha=0) )[[3]]
    lmb2 <- cv2$lambda.min;
    timeFIT2 = system.time( f5 <- glmnet(x=gen[-w,],y=y[-w],lambda = lmb2,alpha = 0) )[[3]]
    ffw5_RidgeReg <- c(predict(f5,gen[w,]))
    timeQ_RidgeReg = timeCV2 + timeFIT2
    #
    timeCV3 = system.time( cv3 <- cv.glmnet(x=gen[-w,],y=y[-w],alpha=0.01) )[[3]]
    lmb3 = cv3$lambda.min;
    timeFIT3 = system.time( f14 <- glmnet(x=gen[-w,],y=y[-w],lambda = lmb3,alpha = 0.01) )[[3]]
    ffw6_ElsNet <- c(predict(f14,gen[w,]))
    timeQ_ElsNet = timeCV3 + timeFIT3
    cat('glmnet\n')
    
    # Boosting
    timeQ_Boosting = system.time( f8 <- gbm.fit(x=gen[-w,],y=y[-w],distribution = "gaussian",verbose = F) )[[3]]
    ff8_Boosting <- predict(f8,gen[w,],100)
    cat('gbm\n')
    
    # REML GBLUP
    timeQ_GREML = system.time( f10 <- mixed.solve(y,K=G) )[[3]] + timeBuildG
    ff110_GREML <- f10$u[w]
    cat('rrBLUP\n')
    
    # PLS
    cat('now PLS\n')
    timeQ_PLS = system.time( f16a <- plsr(y[-w]~gen[-w,],ncomp=5) )[[3]]
    ff16_PLS <- predict(f16a,gen[w,])[,1,5]
    
    # FLM
    timeQ_FLM = system.time( f19b <- FLM(y[-w],gen[-w,]) )[[3]]
    ff19b_FLM <- gen[w,] %*% f19b$b + f19b$mu
    cat('FLM\n')
    
    # EBL
    timeQ_EBL = system.time( f20 <- EBglmnet(x=gen[-w,], y=y[-w],family="gaussian",prior= "lassoNEG",hyperparameters=c(0.5,0.5)) )[[3]]
    B = rep(0,ncol(gen))
    B[f20$fit[,1]] = f20$fit[,3]
    ff20_EmpBL = gen[w,] %*% B + f20$Intercept
    cat('Fast BLASSO1\n')
    
    # Legarra LASSO
    timeQ_ExtBL = system.time( f17d <- vigor(y,gen,'EBL',hyperpara(gen,Mvar=.5,'EBL',.01, f=1),Function = 'tuning',Printinfo = FALSE) )[[3]]
    ff17d_ExtBL = gen[w,] %*% f17d$Beta + f17d$Alpha
    cat('Fast BLASSO2\n')
    
    # Collect prediction output
    NamesMod = grep('^ff.+',ls(),value = T)
    M = data.frame(OBSERVATION = Y[w],mget(NamesMod))
    M = c(cor(M,use='p')[1,-1])
    names(M) = gsub('^.+_','',names(M))
    
    # Collect time output
    NamesMod2 = grep('^timeQ.+',ls(),value = T)
    N = c(mget(NamesMod2))
    names(N) = gsub('^.+_','',names(N))
    
    # Output
    out = list(PA=M,ET=unlist(N))
    return(out)
    
  }
  
  # Running loop
  if(is.null(SET)){
    b = list()
    for(i in 1:cvs){b[[i]] = folds(Seeds[i]); cat('Done with',i,'of',cvs,'\n')}
    file.remove(list = list.files(pattern = '\\.dat'))
    names(b) = paste('CV_',1:length(b),sep='')
  }else{
    b = list()
    for(i in unique(SET)){b[[i]] = folds(Seeds[i],which(SET==i)); cat('Done with',i,'\n')}
    file.remove(list = list.files(pattern = '\\.dat'))
    ww = which(!sapply(b,is.null )); b = b[ww]
    names(b) = paste('CV_',unique(SET),sep='')
  }
  
  # Prepare output
  
  # Elapsed time
  ET = t(sapply(b, function(x) x$ET ))
  ET = apply(ET,2,as.numeric)
  # Predictive ability
  PA = t(sapply(b, function(x) x$PA ))
  PA = apply(PA,2,as.numeric)
  
  # Sort boxplot
  if(box){
    ET = ET[,names(sort(colMeans(ET,na.rm=T),F))]
    PA = PA[,names(sort(colMeans(PA,na.rm=T),T))]
  }else{
    ET = colMeans(ET,na.rm=T)
    PA = colMeans(PA,na.rm=T)
  }
  
  # Return output
  out = list(ET=ET,PA=PA)
  return(out)
}

# Small example
RunExample = FALSE
if(RunExample){
  data(tpod,package = 'NAM')
  H = CV(y,gen,3,1:3)
  par(mfrow=c(2,1))
  boxplot(H$PA,main='Predictive ability',ylab='Cor( Obs, Pred )',las=2)
  boxplot(H$ET,main='Elapsed time',ylab='Seconds',las=2)
}

