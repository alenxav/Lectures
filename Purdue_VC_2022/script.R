# Get some data
tr = function(x) sum(diag(x))
data(met,package='NAM')
Obs$Sp = with(Obs,NAM::SPC(YLD,Block,Row,Col,4,4))
Obs$Block = factor(Obs$Block)
plot(Obs$Sp,Obs$YLD)
Gen = data.matrix(Gen[grep('-05|-15',rownames(Gen)),])
Obs = droplevels(Obs[grep('-05|-15',Obs$ID),])
Gen = Gen[,apply(Gen,2,var)>0.1]
dim(Gen)
dim(Obs)

# Genomic Relationship Matrix
N = apply(Gen,2,function(x)x-mean(x))
A = tcrossprod(N)
SumVarSNPs = mean(diag(A))
A = A/SumVarSNPs
diag(A) = diag(A)+0.01 # Stabilize GRM
A = A/mean(diag(A))
iA = solve(A)

# Design matrices
y = Obs$YLD
X = model.matrix(~Block+Year:Sp-1,Obs)
Z = model.matrix(~ID-1,Obs)

# Constants
n = length(y) # number of obs
q = ncol(Z) # levels of Z
rX = ncol(X) # rank of X

# Starting values for variance components
b = qr.solve(X,y)
vy0 = c(crossprod(y-X%*%b)/(n-rX))
vu = 0.25*vy0
ve = 0.75*vy0
vc = c(vu=vu,ve=ve)
print(vc)


# Variance component matrices
I = diag(n)
R = I*ve
G = A*vu
iG = solve(G)
iR = solve(I*ve)
ZAZ = Z %*% A %*% t(Z)
V = ZAZ*vu + R
iV = solve(V)  

# Mixed model equation matrices
Sigma = Matrix::bdiag(diag(0,rX),iG)
W = cbind(X,Z)
C = t(W) %*% iR %*% W + Sigma
iC = solve(C)
r = t(W) %*% iR %*% y

# Coefficients
g = iC %*% r
b = g[1:rX]
u = g[-c(1:rX)]


C22 = iC[-c(1:rX),-c(1:rX)]
C22 = as.matrix(C22)
S = I - X %*% solve( t(X)%*%X ) %*% t(X)
M = I - X %*% solve( t(X) %*% iV %*%X ) %*% t(X) %*% iV
P = iV - iV %*% X %*% solve( t(X)%*%iV%*%X ) %*% t(X) %*% iV
# Also # P = iV %*% M
# Also # H = iR %*% (I - W %*% iC %*% t(W) %*% iR)

# Variance of u hat
VarUhat = G %*% t(Z) %*% P %*% Z %*% G
VarUhat2 = G - C22
plot(diag(VarUhat),diag(VarUhat2))
# Reliability
rel1 = sqrt(diag(diag(q)-C22 %*% iG))
rel2 = sqrt(diag( t(Z) %*% iV %*% Z %*% G ))
plot(rel1,rel2)

# Other ways to get BLUPs
u1 = g[-c(1:rX)]
u2 = G %*% t(Z) %*% P %*% y
u3 = solve( t(Z)%*%iR%*%Z + iG , t(Z)%*%iR%*%c(y-X%*%b) )
u4 = solve( t(Z)%*%Z + iA*(ve/vu) , t(Z)%*% c(y-X%*%b) )
plot(data.frame(u1,u2,u3,u4))

# Null space projection
yHat = X%*%b + Z%*%u
e = y - yHat
plot(P%*%y*ve,e)

############
# Check VC #
############

# EM
# Var e
Ve = c(crossprod(y,e)) / (n-rX)
# Var U
Vu = c(t(u) %*% iA %*% u + tr(iA%*%C22)) / q
# Var U (faster converging alternative)
Vu2 = c(t(u) %*% iA %*% u) / ( q - tr(iA%*%C22)/vu )
# Check
print(c(Ve=Ve, Vu=Vu, Vu2=Vu2))

# Schaeffer's Pseudo-Expectation
Sy = S %*% y
ZSy = t(Z) %*% Sy
trSZAZ = tr( S %*% ZAZ ) # tr(ZS'Z) if A=I
# Var U
c(u %*% ZSy) / trSZAZ
# Var E
c(t(e) %*% Sy) / (n-rX)

# MIVQUE
lhs = matrix(c(
  tr( P %*% ZAZ %*% P %*% ZAZ ), tr( P %*% ZAZ %*% P %*% I ),
  tr( P %*% I %*% P %*% ZAZ ),  tr( P %*% I %*% P %*% I ) ),2,2)
rhs = c(t(y) %*% P %*% ZAZ %*% P %*% y,
        t(y) %*% P %*%  I %*% P %*% y)
solve(lhs,rhs)

# AI  via V

SecDer1 = matrix(c(
  t(y) %*% P %*% ZAZ %*% P %*% ZAZ %*% P %*% y, t(y) %*% P %*% ZAZ %*% P %*% I %*% P %*% y,
  t(y) %*% P %*% I %*% P %*% ZAZ %*% P %*% y, t(y) %*% P %*% I %*% P %*% I %*% P %*% y
),2,2)
FirDer1 = c( vu = tr(P %*% ZAZ) - t(y) %*% P %*% ZAZ %*% P %*% y ,
             ve = tr(P %*% I) - t(y) %*% P %*% I %*% P %*% y )
vc - solve(SecDer1,FirDer1)

# AI via C

#B = cbind( vu = ZAZ %*% P %*% y, ve = I %*% P %*% y)
B = cbind( vu = Z %*% u / vu, ve = e / ve )

MB = chol(rbind(cbind(as.matrix(C),       t(W) %*% iR %*% B),
                cbind( t(B) %*% iR %*% W, t(B) %*% iR %*% B) ))

LB = MB[(ncol(MB)-1):ncol(MB),(ncol(MB)-1):ncol(MB)]
SecDer2 = crossprod(LB)
FirDer2 = c( vu = ( q/vu - (t(u)%*%iA%*%u)/(vu^2) - tr(iA%*%C22)/(vu^2) ),
             ve = ( (n-rX)/ve - (1/ve)*(q-tr(iA%*%C22)/(vu))) - crossprod(e)/(ve^2) )
vc - solve(SecDer2,FirDer2)

# Gibbs sampler

# Priors
df0 = 5
Su = vu
Se = ve

# Samples of VC
(t(u) %*% iA %*% u + Su*df0 ) / rchisq(1,q+df0)
(t(e) %*% e + Se*df0 ) / rchisq(1,n+df0)

# Sample coefficients
g_samp = sapply(1:(rX+q), function(j) rnorm(1,mean=c(r[j]-C[j,-j]%*%g[-j])/C[j,j],sd=sqrt(1/C[j,j])))
plot(g,g_samp)


# Marker effects
vb = vu/SumVarSNPs
beta = vb * c( t(Z1) %*% iG %*% u)
plot(beta)

### GAUSS-SEIDEL

# Get some data
data(tpod, package='NAM')
y = y-mean(y) # centralize phenotype
X = NAM::CNT(gen) # centralize markers
b = rep(0,ncol(X)) # coefficient starting values
xx = apply(X,2,crossprod) # pre-compute X'X
lambda = mean(xx) # lambda for h2=0.5
e = y # starting value of e

# Gauss-Seidel
for(j in sample(1:ncol(X))){
  b_old = b[j]
  b[j] = (c(X[,j]%*%e) +xx[j]*b_old)/(xx[j]+lambda)
  e = e - X[,j]*(b[j]-b_old)
}

# Check convergence over 10 iterations
for(i in 1:10){
  # Store current version of beta
  b_vec = b
  # Gauss-Seidel
  for(j in sample(1:ncol(X))){
    b_old = b[j]
    b[j] = (c(X[,j]%*%e) +xx[j]*b_old)/(xx[j]+lambda)
    e = e - X[,j]*(b[j]-b_old)}
  # Print log convergence
  cat(round(log10(crossprod(b_vec-b)[1,1]),2),' ')
}


plot(b,xlab='SNP',ylab='Marker effect')

