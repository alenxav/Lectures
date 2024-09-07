load("SAMPLE_G2F4GXE.RData")
require(bWGR)

# Prep matrix for MegaLMM
YY = data.frame(Y_ES)[rownames(gen2),]
rownames(YY) = rownames(gen2)
YY = as.matrix(YY)

# PC matrix
ZZp = EigenGRM(gen2)
system.time(X<-K2X(ZZp))[3]
v = apply(X,2,var)
X = X[,which(v>1e-7)]
rownames(X) = rownames(gen2)
i = intersect(rownames(Y_ES),rownames(X))
j = intersect(rownames(Y_PS),rownames(X))

# Validation function
validation = function(hat,avg=TRUE){
  colnames(hat) = colnames(Y_ES)
  st1 = substr(colnames(hat),1,4)
  st2 = substr(colnames(Y_PS),1,4)
  prd = t(apply(hat,1,function(x) tapply(x,st1,mean))[st2,])
  x = data.frame(
    ByState = diag(cor(prd,Y_PS,use='p')),
    AllxAll =  colMeans(cor(hat,Y_PS,use='p')),
    Average =  cor(rowMeans(hat),Y_PS,use='p')[1,])
  xx = apply(x,2,function(xx){
    xx = round(c(Mu=mean(xx),Std=sd(xx)),2)
    paste0(xx[1],' (',xx[2],')')}) 
  if(avg) return(xx) else return(x)
}

# UV across
y2 = rowMeans(apply(Y_ES,2,function(x) x-mean(x,na.rm=T)),na.rm = T)
y2 = y2[!is.na(y2)]
ii = intersect(names(y2),rownames(X))
uv0 = emDE(scale(y2[ii]),X[ii,])$b
mean(cor(X[j,] %*% uv0,Y_PS ,use='p'))

# MV on Residuals
suv0 = c(X[i,] %*% uv0)
R_ES = apply(Y_ES[i,],2,scale)-suv0
sem1 = ZSEMF(R_ES,X[i,])$b
hat = X[j,] %*% (sem1+uv0)
validation(hat)
