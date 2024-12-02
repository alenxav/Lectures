require(bWGR)

# Simulate data
Z = SimZ()
GC = SimGC()
S = SimY(Z=Z, GC=GC)
Y = S$Y
tbv = S$tbv

# Fit models
system.time(fit_Univariate <- Z %*% FUVBETA(Y,Z))[3]
system.time(fit_THGS <- MRR3F(Y,Z,TH=TRUE)$hat)[3]
system.time(fit_MegaSEM <- SEM(Y,Z)$hat)[3]

# Check accuracy
sort(sapply(grep('fit',ls(),value=T),function(x) mean(diag(cor(get(x),tbv)))))
