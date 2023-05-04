
# Check if you have all the packages, install otherwise

Pkgs = c("NAM","bWGR","ranger","readr","magrittr") 
MyInstalledPkgs = installed.packages()[, 1]                                                                                                                   
check = function(mypkg){ if(!is.element(mypkg, MyInstalledPkgs)) {install.packages(mypkg, repos = "https://cloud.r-project.org")}}
sapply(Pkgs, check)
rm(Pkgs,MyInstalledPkgs,check)

###############################
# PREPARE GENOMIC INFORMATION #
###############################

# Read in genomic info
require(readr)
gen = read_tsv('Training_Data/5_Genotype_Data_All_Years.vcf')
genID = gen$ID
gen2 = gen[,-c(1:9)]
rm(gen)
gc()

# CHUNK ADDED TO TEST SCRIPT TO SPEED UP SCRIPT VALIDATION
set.seed(123)
sample_geno = sort(sample(nrow(gen2),35000))
gen2 = gen2[sample_geno,]
genID = genID[sample_geno]
gc()

# Formate genomic info
for(i in 1:ncol(gen2)){
  cat(i,'\n')
  x = gen2[,i][[1]]
  x[x=='0/0']=0
  x[x=='0/1'|x=='1/0']=1
  x[x=='1/1']=2
  x = as.numeric(x)
  gen2[,i] = x}
rm(i,x)
gc()

require(magrittr)
gen2 %<>% as.data.frame()
rownames(gen2) = genID
gc()

for(i in 1:ncol(gen2)){ cat(i,'\n'); gen2[,i] = as.integer(gen2[,i])}
rm(i)
gc()

# Format correctly
gen2 %<>% as.matrix()
gc()
gen2 %<>% t()
gc()
genID = rownames(gen2)

# remove odd behaving SNPs
m = apply(gen2,2,function(x) max(x,na.rm = T))
m = which(m==1|m==2)
gen2 = gen2[,m]
rm(m)
gc()

# drop markers with too much missing
m = apply(gen2,2,function(x) mean(!is.na(x)) )
gen2 = gen2[,m>0.9]
gc()

# drop markers based on MAF
PS_dta = read.csv('Testing_Data/1_Submission_Template_2022.csv')
pt = unique(PS_dta$Hybrid)
w = which(rownames(gen2)%in%pt)
m = colMeans(gen2[w,],na.rm = T) * 0.5 # PS MAF
gen2 = gen2[,m>0.1]
gc()
m = colMeans(gen2[-w,],na.rm = T) * 0.5  # ES MAF
gen2 = gen2[,m>0.1]
gc()

# imputation
m = c(table(gsub('_.+','',colnames(gen2))))
gen2 = bWGR::markov(gen2,m)
gc()

# drop redundant markers (full LD or close to it)
set.seed(123)
m = matrix(rnorm(3*nrow(gen2)),ncol=3)
ev = t(gen2) %*% m
unique_ev = unique(round(ev))
gen2 = gen2[,rownames(unique_ev)]
rm(m,ev,unique_ev)

# multiple rounds of dropping neighbor SNPs based on LD
gen2 = NAM::snpQC(gen2,psy=0.999)          # 1st round
gen2 = NAM::snpQC(gen2,psy=0.999, MAF = 0) # 2nd round
gen2 = NAM::snpQC(gen2,psy=0.998, MAF = 0) # 3rd round
gen2 = NAM::snpQC(gen2,psy=0.998, MAF = 0) # 4th round
gen2 = NAM::snpQC(gen2,psy=0.997, MAF = 0) # 5th round
gen2 = NAM::snpQC(gen2,psy=0.997, MAF = 0) # 6th round
gen2 = apply(gen2,2,as.integer)
rownames(gen2) = genID
gc()
save(gen2, file='gen.RData')

# Get kernel
K = bWGR::EigenARC(gen2,T,4)
rownames(K) = colnames(K) = rownames(gen2)
save(K, file = 'ArcKernel.RData')
rm(gen2); gc()

# Reparametrize
require(bWGR)
K2X = function(K){
  E = EigenEVD(K,cores=4)
  w = which(E$D>0.1)
  X = E$U[,w] %*% diag(sqrt(E$D[w]))
  rownames(X) = rownames(K)
  return(X)}
X = K2X(K)
save(X,file='K2X.RData')

# ----------------------------------------------------------------------------------------------------------------------------------------------------

######################################
# Catch outliers, prepare parameters #
######################################

### Read other files
ES_dta = read.csv('Training_Data/1_Training_Trait_Data_2014_2021.csv')
ES_meta = read.csv('Training_Data/2_Training_Meta_Data_2014_2021.csv')
ES_soil = read.csv('Training_Data/3_Training_Soil_Data_2015_2021.csv')
ES_weat = read.csv('Training_Data/4_Training_Weather_Data_2014_2021.csv')
ES_EC = read.csv('Training_Data/6_Training_EC_Data_2014_2021.csv',row.names = 1)
PS_dta = read.csv('Testing_Data/1_Submission_Template_2022.csv')
PS_meta = read.csv('Testing_Data/2_Testing_Meta_Data_2022.csv')
PS_soil = read.csv('Testing_Data/3_Testing_Soil_Data_2022.csv')
PS_weat = read.csv('Testing_Data/4_Testing_Weather_Data_2022.csv')
PS_EC = read.csv('Testing_Data/6_Testing_EC_Data_2022.csv',row.names = 1)

# Add meta covariates to ES and PS
rownames(ES_meta) = ES_meta$Env
rownames(PS_meta) = PS_meta$Env
ES_meta$State = gsub('H.+|S1.+','',rownames(ES_meta))
PS_meta$State = gsub('H.+|S1.+','',rownames(PS_meta))
ES_meta$State[ES_meta$State=='O'] = 'OH'
PS_meta$State[PS_meta$State=='O'] = 'OH'
unique(PS_meta$State) %in% unique(ES_meta$State)
ES_meta$Irr = grepl('Irr',ES_meta$Treatment)
PS_meta$Irr = PS_meta$Irrigated == 'yes'
ES_dta$Irr = ES_meta[ES_dta$Env,'Irr']
PS_dta$Irr = PS_meta[PS_dta$Env,'Irr']
ES_dta$Trt = ES_meta[ES_dta$Env,'Treatment']
PS_dta$Trt = PS_meta[PS_dta$Env,'Treatment']
ES_dta$PC = ES_meta[ES_dta$Env,'Previous_Crop']
PS_dta$PC = PS_meta[PS_dta$Env,'Previous_Crop']
ES_dta$State = ES_meta[ES_dta$Env,'State']
PS_dta$State = PS_meta[PS_dta$Env,'State']

# Pooling levels
PS_dta$PC[PS_dta$PC%in%c('peanut','soybean','soybeans with fall cereal rye cover')] = 'legume'
ES_dta$Trt[ES_dta$Trt%in%c("Late planting")] = 'Late'
ES_dta$Trt[ES_dta$Trt%in%c('Dry Land','Dryland','Dryland optimal')] = 'Dry'
ES_dta$Trt[ES_dta$Trt%in%c('Standard - Irrigated Optimal','Early Planting','','Irrigated')] = 'Standard'
ES_dta$Trt[ES_dta$Trt%in%c('Late Planted Irrigated','Late Planting','Late Stressed')] = 'Late'
ES_dta$PC[ES_dta$PC%in%c('peanut','soybean',"soybean/pumpkin","wheat/soybean","wheat and Double Crop soybean",
                         "wheat/double crop soybean","Lima beans followed by rye cover crop","Small Grains and Double Crop soybean")] = 'legume'
ES_dta$PC[ES_dta$PC%in%c("Winter wheat","2019/20 wheat")] = 'wheat'
ES_dta$PC[ES_dta$PC%in%c("sugar beet")] = ''
ES_dta$PC[ES_dta$PC=="Fallow most of 2014 winter planted in fall of 2014 then sprayed with Glystar 24 floz/a on 5/3/15  and killed spring of 2015 spray "]=''
unique(PS_dta$Trt) %in% unique(PS_dta$Trt)

# Catch outliers
ES_dta$ENY = EnvNormYield = (ES_dta$Yield_Mg_ha - tapply(ES_dta$Yield_Mg_ha,ES_dta$Env,mean,na.rm=T)[ES_dta$Env])/tapply(ES_dta$Yield_Mg_ha,ES_dta$Env,sd,na.rm=T)[ES_dta$Env]
ES_dta$YLD = ifelse(abs(ES_dta$ENY)<3,ES_dta$Yield_Mg_ha,NA)
save(ES_dta,PS_dta,file='ES_dta_QC.RData')

# ----------------------------------------------------------------------------------------------------------------------------------------------------

################################
# Modeling Environmental Means #
################################

load('ES_dta_QC.RData')

# Estimation set (training set)
EMD = unique(ES_dta[,c('Env','Irr','Trt','PC','State')])
rownames(EMD) = EMD$Env
EMD = EMD[,-1]
EMD$YLD = tapply(ES_dta$YLD, ES_dta$Env, mean,na.rm=T)[rownames(EMD)]
# Prediction set (test set)
EMDp = unique(PS_dta[,c('Env','Irr','Trt','PC','State')])
rownames(EMDp) = EMDp$Env
EMDp = EMDp[,-1]
EMDp$YLD = NA
# EC (environmental covariates)
X_EC = rbind(ES_EC,PS_EC)
EMD2 = rbind(EMD,EMDp)[rownames(X_EC),]

# fit some model (RF and linear)
require(ranger)
require(bWGR)
EMD3X = cbind(EMD2,X_EC)
EMD3 = EMD3X[!is.na(EMD3X$YLD),]
fit_emd = ranger(YLD~.,EMD3,1500,500,always.split.variables = colnames(EMD)[1:4])

# Reparametrize W
W2X = function(Z,AK=FALSE,Phi=1){
  if(AK){ K = EigenARC(Z)           }else{
          K = EigenGAU(Z,phi = Phi) }
  E = eigen(K,symmetric = T)
  E$values[E$values<0.0001] = 0
  X = E$vectors %*% diag(sqrt(E$values))
  X = X[,which(E$values>0)]
  rownames(X) = rownames(W)
  return(X)}
W = model.matrix(~.-1,EMD3X[,-5])
vW = apply(W,2,var,na.rm=T)
W3 = W2X(W[,which(vW>0)],Phi=2)
W3ES = W3[rownames(ES_EC),]
yr = gsub('^.+_','',rownames(W3ES))
lc = substr(rownames(W3ES),1,2)
fit_emd_ba = emBA(EMD3$YLD,W3ES)

# Simple linear model (least squares)
X0 = model.matrix(~.-1,EMD2[,-5])
X = X0[rownames(ES_EC),]
X1 = X0[rownames(PS_EC),]
y = EMD2[rownames(X),'YLD']
fit_lm = emML(y,X)

# Average hats, minimize bias
H1 = predict(fit_emd,EMD3X)$predictions
H2 = fit_emd_ba$mu + W3 %*% fit_emd_ba$b
H3 = fit_lm$mu + X0 %*% fit_lm$b
Avg = c(H1+H2+H3)/3
Y = EMD3X$YLD
fit_avg = lm(Y~Avg)
Avg2 = predict(fit_avg,data.frame(Avg=c(Avg)))
DF = data.frame(Observed=Y,RndForst=H1,BayesA_EC=H2,LinearModel=H3,AvgModel=Avg2)

# Filling missing predictions with averages of previous years
if(any(!unique(PS_dta$Env)%in%rownames(DF))){
  missingEnv = unique(PS_dta$Env)[(!unique(PS_dta$Env)%in%rownames(DF))]
  codes = substr(missingEnv,1,4)
  obs_codes = substr(rownames(DF),1,4)
  imp = t(sapply(codes, function(x) colMeans(DF[obs_codes==x,])))
  rownames(imp) =  missingEnv
  DF = rbind(DF,imp)}

write.csv(DF,'EM_models.csv')

# ----------------------------------------------------------------------------------------------------------------------------------------------------

#############################
# Modeling Environmental SD #
#############################

#load("image_with_env_mean_models.RData")
EMD$YLDsd = tapply(ES_dta$YLD, ES_dta$Env, sd,na.rm=T)[rownames(EMD)]
EMDp$YLDsd = NA
EMD2SD = rbind(EMD,EMDp)[rownames(X_EC),]
EMD3XSD = cbind(EMD2SD,X_EC)
EMD3SD = EMD3XSD[!is.na(EMD3XSD$YLD),]
EMD3XSD = EMD3XSD[,which(names(EMD3XSD)!='YLD')]
EMD3SD = EMD3SD[,which(names(EMD3SD)!='YLD')]

require(ranger)
fit_emd_sd = ranger(YLDsd~.,EMD3SD,1500,500)
Hsd = predict(fit_emd_sd,EMD3XSD)$predictions
names(Hsd) = rownames(EMD3XSD)

# Impute missing
PS_dta = read.csv('Testing_Data/1_Submission_Template_2022.csv')

if(any(!unique(PS_dta$Env)%in%names(Hsd))){
  missingEnv = unique(PS_dta$Env)[(!unique(PS_dta$Env)%in%names(Hsd))]
  codes = substr(missingEnv,1,4)
  obs_codes = substr(rownames(DF),1,4)
  imp = sapply(codes, function(x) mean(Hsd[obs_codes==x],na.rm=T))
  names(imp) =  missingEnv
  Hsd = c(Hsd,imp)}
write.csv(Hsd,file='ESD.csv')

# ----------------------------------------------------------------------------------------------------------------------------------------------------

######################
# Compute Accuracies #
######################

load('K2X.RData')
load('ES_dta_QC.RData')

# identify individuals from ES and PS
es = tapply(ES_dta$Hybrid, ES_dta$Env, function(x) unique(x[x%in%rownames(X)])    )
ps = tapply(PS_dta$Hybrid, PS_dta$Env, function(x) unique(x[x%in%rownames(X)])    )

# Get accuracy between every combination of training and test set
A = matrix(NA,length(es),length(ps),dimnames = list(names(es),names(ps)))
require(bWGR)
for(i in names(es)){
  for(j in names(ps)){
    if(is.na(A[i,j])){
      cat(i,length(es[[i]]),'|',j,length(ps[[j]]),'\n')
      A[i,j] = mean(EigenAcc(X[es[[i]],],X[ps[[j]],]))
      cat(A[i,j],'\n')
    }
  }  
}

save(A,file='AccESPS.RData')

# ----------------------------------------------------------------------------------------------------------------------------------------------------

##################
# Modeling GEBVs #
##################

require(bWGR)
load('ES_dta_QC.RData')
if(!'K'%in%ls()) load('ArcKernel.RData')

AdjYld = list()
envs = unique(ES_dta$Env)

# Spatial adjustment function
SpAdj = function(env){
  cat('Location: ',env,'\n')
  dta = subset(ES_dta,Env == env)
  cat('Obs ',nrow(dta),'\n')
  dta = droplevels(dta[dta$Hybrid%in%rownames(K),])
  id = unique(dta$Hybrid)
  cat('Hyb ',length(id),'\n')
  trt = c('YLD','Stand_Count_plants','Pollen_DAP_days','Silk_DAP_days','Plant_Height_cm','Ear_Height_cm')
  trt = trt[sapply(dta[trt],function(x) mean(!is.na(x)) > 0.5 & var(x,na.rm=T)>0 )]
  cat('Trt ',length(trt),'\n')
  if(!anyNA(dta$Pass)){
    Q = NAM::covar(dta[,c('Block','Range','Pass')],rho=3.5,type=1.75,dist=1.2)
    G = K[dta$Hybrid,dta$Hybrid]
    Y = data.matrix(dta[,trt])
    fit = tryCatch(mkr2X(Y,G,Q),error = function(e) NULL)
    if(is.null(fit)){
      cat('No spatial, convergence problem\n');
      y = c(dta[,'YLD'])
    }else{
      y = Y[,1]-fit$b2[,1];
      # Checks
      if(all(is.na(y))){
        cat('Algorithm broke, use pheno\n');
        y = c(dta[,'YLD'])
      }else if(mean(y,na.rm=T)<0){
        cat('Negative mean (odd results, use pheno)\n');
        y = c(dta[,'YLD'])
      }else{
        tmp = round(fit$Vb2[1,1]/var(Y[,1],na.rm=T),2)
        cat('Sp var removed',tmp,'\n')
        if(tmp<0|tmp>1){y = c(dta[,'YLD']);
        cat('Odd results, use pheno \n')
      }else{
        cat('Spatial adjustment successful\n')
        }
      }
    } 
  }else{ cat('No spatial information\n'); y = c(dta[,'YLD']) }
  yy = tapply(y,dta$Hybrid,mean,na.rm=T)
  yy = yy[rownames(K)]
  cat('Done\n\n')
  return(yy)
}

for(i in envs){
  if(!i %in% names(AdjYld)){
    cat(length(AdjYld)+1,'\n')
    AdjYld[[i]] = SpAdj(i)
    save(AdjYld,file = 'AdjYld.RData')
  }
}

# Fit MV with spatially adjusted phenotypes
require(bWGR)
if(!'AdjYld'%in%ls()) load('AdjYld.RData')
if(!'X'%in%ls()) load('K2X.RData')

Y = sapply(AdjYld,c)
rownames(Y) = rownames(X)

if('ES_MVGBLUP1.RData'%in%dir()){
  load('ES_MVGBLUP1.RData')
}else{
  fit_mv = MRR3(Y,X,InnerGS=T,NoInv=T,df0=10,cores=4,tol=1e-04)
  fit_mv$Its
  colnames(fit_mv$GC) = rownames(fit_mv$GC) = 
    colnames(fit_mv$hat) = colnames(fit_mv$b) = 
    names(fit_mv$h2) = colnames(Y)
  if(all(!is.na(fit_mv$h2))) save(fit_mv,Y,file = 'ES_MVGBLUP1.RData')
}

# ----------------------------------------------------------------------------------------------------------------------------------------------------

###############################
# Final predictions - Model A #
###############################


# Get submission template
PS_dta = read.csv('Testing_Data/1_Submission_Template_2022.csv')

# Get pre-calculated environmental means and std dev
EPRED = read.csv('EM_models.csv',row.names = 1)
Hsd = rowMeans(read.csv('ESD.csv',row.names = 1))
# Compute genomic values
if(!'X'%in%ls()) load('K2X.RData')
if(!'fit_mv'%in%ls()) load('ES_MVGBLUP1.RData')
H = X %*% fit_mv$b
H = apply(H,2,scale)

load('AccESPS.RData') # load accuracies between ES and PS
A = A[colnames(H),]
ES_state = substr(rownames(A),1,2)
PS_state = substr(colnames(A),1,2)
ES_code = substr(rownames(A),1,4)
PS_code = substr(colnames(A),1,4)

# Selection index
SI  = A
for(i in 1:ncol(SI)){
  # Index as a function of code, state and relationship ESPS
  index = (ES_code == PS_code[i])*A[,i]*2 + (ES_state == PS_state[i])*A[,i]*2 + A[,i]*0.1
  index = index/sum(index)
  SI[,i] = index}

GPRED = H %*% SI
rownames(GPRED) = rownames(X)

for(i in colnames(GPRED))  GPRED[,i] = GPRED[,i]*Hsd[i]
write.csv(GPRED,'GM_models_PredVar.csv')

GPRED[1:5,1:5]
EPRED[1:5,1:5]

hat2022 = apply(PS_dta,1,function(x) GPRED[x[2],x[1]]+EPRED[x[1],'AvgModel']  )
PS_dta$Yield_Mg_ha = hat2022
write.csv(PS_dta,'MODEL_A.csv', row.names = F, quote = F)

# ----------------------------------------------------------------------------------------------------------------------------------------------------

##############################
# Second model: Simple GBLUP #
##############################

# LOAD DATA WITH NO QC
ES_dta = read.csv('Training_Data/1_Training_Trait_Data_2014_2021.csv')
ES_meta = read.csv('Training_Data/2_Training_Meta_Data_2014_2021.csv')
PS_dta = read.csv('Testing_Data/1_Submission_Template_2022.csv')
PS_meta = read.csv('Testing_Data/2_Testing_Meta_Data_2022.csv')

# Add meta covariates to ES and PS
rownames(ES_meta) = ES_meta$Env
rownames(PS_meta) = PS_meta$Env
ES_meta$State = gsub('H.+|S1.+','',rownames(ES_meta))
PS_meta$State = gsub('H.+|S1.+','',rownames(PS_meta))
ES_meta$State[ES_meta$State=='O'] = 'OH'
PS_meta$State[PS_meta$State=='O'] = 'OH'
unique(PS_meta$State) %in% unique(ES_meta$State)
ES_meta$Irr = grepl('Irr',ES_meta$Treatment)
PS_meta$Irr = PS_meta$Irrigated == 'yes'
ES_dta$Irr = ES_meta[ES_dta$Env,'Irr']
PS_dta$Irr = PS_meta[PS_dta$Env,'Irr']
ES_dta$Trt = ES_meta[ES_dta$Env,'Treatment']
PS_dta$Trt = PS_meta[PS_dta$Env,'Treatment']
ES_dta$PC = ES_meta[ES_dta$Env,'Previous_Crop']
PS_dta$PC = PS_meta[PS_dta$Env,'Previous_Crop']
ES_dta$State = ES_meta[ES_dta$Env,'State']
PS_dta$State = PS_meta[PS_dta$Env,'State']
ES_dta$Station = substr(ES_dta$Env,1,4) # added on 12/26
PS_dta$Station = substr(PS_dta$Env,1,4) # added on 12/26
unique(PS_dta$Station) %in% unique(ES_dta$Station)
# Pooling levels
PS_dta$PC[PS_dta$PC%in%c('peanut','soybean','soybeans with fall cereal rye cover')] = 'legume'
ES_dta$Trt[ES_dta$Trt%in%c("Late planting")] = 'Late'
PS_dta$Trt[PS_dta$Trt%in%c("Late planting")] = 'Late' # fixed on 12/25
ES_dta$Trt[ES_dta$Trt%in%c('Dry Land','Dryland','Dryland optimal')] = 'Dry'
ES_dta$Trt[ES_dta$Trt%in%c('Standard - Irrigated Optimal','Early Planting','','Irrigated')] = 'Standard'
ES_dta$Trt[ES_dta$Trt%in%c('Late Planted Irrigated','Late Planting','Late Stressed')] = 'Late'
ES_dta$PC[ES_dta$PC%in%c('peanut','soybean',"soybean/pumpkin","wheat/soybean","wheat and Double Crop soybean",
                         "wheat/double crop soybean","Lima beans followed by rye cover crop","Small Grains and Double Crop soybean")] = 'legume'
ES_dta$PC[ES_dta$PC%in%c("Winter wheat","2019/20 wheat")] = 'wheat'
ES_dta$PC[ES_dta$PC%in%c("sugar beet")] = ''
PS_dta$PC[PS_dta$PC%in%c("sorghum")] = ''# fixed on 12/25
ES_dta$PC[ES_dta$PC=="Fallow most of 2014 winter planted in fall of 2014 then sprayed with Glystar 24 floz/a on 5/3/15  and killed spring of 2015 spray "]=''
unique(PS_dta$Trt) %in% unique(ES_dta$Trt)
ES_dta$PC[ES_dta$PC%in%c("")] = 'UKNOWN'
PS_dta$PC[PS_dta$PC%in%c("")] = 'UKNOWN'

tapply(ES_dta$Env, ES_dta$Irr, function(x) length(unique(x)))
tapply(PS_dta$Env, PS_dta$Irr, function(x) length(unique(x)))

require(bWGR)
if(!'X'%in%ls()) load('K2X.RData')
Both = rbind(cbind(PS_dta[,c('Env','Hybrid','Irr','Trt','PC','State','Station')], Yield_Mg_ha=NA),
             ES_dta[,c('Env','Hybrid','Irr','Trt','PC','State','Station','Yield_Mg_ha')])

Both$Irr = as.numeric(Both$Irr)
Both$Hybrid = factor(Both$Hybrid)
Both$Trt = factor(Both$Trt)
Both$PC = factor(Both$PC)
Both$State = factor(Both$State)
Both$Station = factor(Both$Station)

# Fit model
fit = mixed(y = Yield_Mg_ha,
            random = ~Hybrid,
            fixed = ~Irr+Trt+PC+State+Station,
            data = Both,
            X = list(Hybrid=X))

S42 = PS_dta
S42$Yield_Mg_ha = fit$Coefficients$Intercept +
  ifelse(S42$Irr,fit$Coefficients$Irr,0) +
  fit$Coefficients$Trt[as.character(S42$Trt)] +
  fit$Coefficients$PC[as.character(S42$PC)] +
  fit$Coefficients$State[as.character(S42$State)] +
  fit$Coefficients$Station[as.character(S42$Station)] +
  c(X[as.character(S42$Hybrid),] %*% fit$Structure$Hybrid) 
S42 = S42[,1:3]
anyNA(S42)

write.csv(S42,'MODEL_B.csv', row.names = F, quote = F)

# ----------------------------------------------------------------------------------------------------------------------------------------------------

#################
# MODEL AVERAGE #
#################


MA = read.csv('MODEL_A.csv')
MB = read.csv('MODEL_B.csv')

FINAL = 0.5 * (MA$Yield_Mg_ha + MB$Yield_Mg_ha)
FF = MA
FF$Yield_Mg_ha = FINAL

write.csv(FF,'CLACsubmission5.csv', row.names = F, quote = F)

