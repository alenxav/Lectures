
# Packages required: lme4, BGlR, asreml (or rrBLUP), eMM2 (or bWGR)

# Load data
load("simulated_set.RData")

# Simulator
SingleStep_simulator = function(tbv, Gen, n = NULL, locs = NULL, SEED = 10){
  
  # Check if you have asreml and eMM2
  AS = 'asreml'%in%rownames(installed.packages())
  EM = 'eMM2'%in%rownames(installed.packages())
  
  # Random values if n, h2 and number of locs are not provided
  set.seed(SEED)
  if(is.null(n)) n = sample(c(100,250,500,1000,1500,2000),1)
  if(is.null(locs)) locs = rpois(1,10)
  h2 = runif(locs,0.25,0.75)
  
  # Produce dataset
  cat('\n','number of observations =',n,'\n','h2 of each loc =',round(h2,2),'\n','seed =',SEED,'\n\n')
  ID = sample(names(tbv),n,replace = T)
  ID0 = factor(ID,levels = names(tbv))
  ENV = factor(1:locs)[sample(locs,n,replace = T)]
  EnvEff = rnorm(locs,0,10)
  Res = rnorm(n,((1-h2[as.numeric(ENV)])/h2[as.numeric(ENV)])*sd(tbv))
  Y = 100 + tbv[ID] + EnvEff[as.numeric(ENV)] + Res
  DTA = data.frame(Y,ID,ENV)
  
  # FLMSS
  if(EM){
    require(eMM2)
    at = system.time(a <- FLMSS(y=Y,random=~ID,fixed=~ENV,data=DTA,M=list(ID=Gen)))[[3]]
    gav = Gen%*%a$Mrk$ID
  }else{
    require(bWGR)
    at = system.time(a <- mixed(y=Y,random=~ID,fixed=~ENV,data=DTA,X=list(ID=Gen),alg=emDE,maxit=5))[[3]]
    gav = Gen%*%a$Structure$ID
  }
  
  # GREML
  timeK = system.time(ZZ <- tcrossprod(apply(Gen,2,function(x) x-mean(x) )))[[3]]
  timeK = system.time(K <- ZZ/mean(diag(ZZ)))[[3]] + timeK
  diag(K) = diag(K)+1e-8
  timeK = system.time(iK <- solve(K))[[3]] + timeK
  rownames(iK)=colnames(iK)=rownames(Gen)
  
  if(AS){
    require(asreml)
    bt = system.time(b <- asreml(Y~ENV,~ped(ID),data=DTA,ginverse=list(ID=iK)))[[3]] + timeK
    gbv = b$coefficients$random
  }else{
    require(rrBLUP)
    Z = model.matrix(~ID0-1); colnames(Z)=gsub('ID','',colnames(Z))
    z = colnames(Z)
    bt = system.time(b <- mixed.solve(y=Y,X=model.matrix(~ENV),Z=Z,K=K))[[3]]
    gbv = b$u
  }
  
  # TWOSTEP
  require(lme4) # For Two-steps (Step1)
  require(BGLR) # For Two-steps (Step2)
  ct = system.time(c1 <- lmer(Y~ID+(1|ENV)-1,data=DTA))[[3]] +
       system.time(wts <- 1/summary(c1)$coefficients[,2] )[[3]] +
       system.time(c2 <- BGLR(c1@beta,ETA=list(g=list(X=Gen[levels(c1@frame$ID),],model='BL')),verbose=F,weights=wts))[[3]]
  gcv = Gen%*%c2$ETA$g$b
  
  # Collect accuracies
  accA = cor(tbv,gav,use='p')
  accB = cor(tbv,gbv,use='p')
  accC = cor(tbv,gcv,use='p')
  cat('\n Accuracy:',paste(c('FLMSS','GBLUP','TWO-STEPS'),round(c(accA,accB,accC),3)),'\n\n')
  
  # Output
  ET = c(FLMSS=at, GBLUP=bt, TWOSTEPS=ct)
  PA = c(FLMSS=accA, GBLUP=accB, TWOSTEPS=accC)
  
  # Return
  out=list(ET=ET,PA=PA)
  return(out)}

# Small example
RunExample = FALSE
if(RunExample){
  TEST = SingleStep_simulator(TBV100,Gen)
  par(mfrow=c(2,1))
  barplot(TEST$PA,ylab='Cor( Pred, TBV )',main='Accuracy')
  barplot(TEST$ET,ylab='Seconds',main='Elapsed time')
}

