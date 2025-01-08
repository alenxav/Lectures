# Code based on the 2022 model

# Read inputs
Z = read.csv('Maize_GxE_Competition_Data_2024/Training_data/5_Genotype_Data_All_2014_2025_Hybrids_numerical.txt', sep='\t', skip=1, row.names=1) |> as.matrix()
dta = read.csv('Maize_GxE_Competition_Data_2024/Training_data/1_Training_Trait_Data_2014_2023.csv')
PSet = read.csv('Maize_GxE_Competition_Data_2024/Testing_data/1_Submission_Template_2024.csv')
require(bWGR) 

# Minor QC
dta$Loc = factor(substr(dta$Env,1,2))
PSet$Loc = substr(PSet$Env,1,2)
dta = subset(dta, Loc %in% unique(PSet$Loc))
SnpVarPSet = apply(Z[unique(PSet$Hybrid),],2,var,na.rm=T)
SnpVarESet = apply(Z[intersect(unique(dta$Hybrid),rownames(Z)),],2,var,na.rm=T)
Z = IMP(Z[,which(SnpVarPSet>0.01 & SnpVarESet>0.01)])
X = K2X(EigenGAU(Z))
rownames(X) = rownames(Z)
colnames(X) = paste0('pc',1:ncol(X))
loc_mu = with(dta, tapply(Yield_Mg_ha,Env,mean,na.rm=T))
dta$cnt = with(dta, Yield_Mg_ha - loc_mu[Env]) 
Y = reshape2::acast(dta, Hybrid~Env, mean, value.var = 'cnt')
Y = Y[,which(colSums(!is.na(Y))>100)]
overlap = intersect(rownames(Y),rownames(X))
dta = subset(dta, Env %in% colnames(Y) & Hybrid %in% overlap)
Y = Y[overlap,]

# Fit single-step model
fit1 = mm( y=Yield_Mg_ha, random=~Env+Hybrid, data=dta,  M=list(Hybrid=X[overlap,]))

# Fit multivariate
fit2 = ZSEMF(Y, X[overlap,], npc = -1)

# Model average and Prediction
B1 = fit1$Mrk$Hybrid 
B2 = t(apply(fit2$b,1,function(x)tapply(x,substr(colnames(Y),1,2),mean) ))
hat = X %*% (B1+B2)

# Validation set
PSet$Loc = substr(PSet$Env,1,2)
colnames(hat) = gsub('Loc','',colnames(hat))
hat = apply(hat,2,function(x) x-mean(x))
mu = tapply(loc_mu,substr(names(loc_mu),1,2),mean)
for(i in 1:nrow(PSet)){  PSet$Yield_Mg_ha[i] = hat[ PSet$Hybrid[i], PSet$Loc[i]] + mu[PSet$Loc[i]] }
PSet = PSet[,1:3]

# Write out CSV
write.csv(PSet, 'Prediction_Submission_Template.csv', row.names = F)

# Check GxE
gc = fit2$GC
rownames(gc) = colnames(gc) = colnames(Y)
heatmap(gc)
