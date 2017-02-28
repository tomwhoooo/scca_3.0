\name{scca-package}
\alias{scca-package}
\alias{scca}
\docType{package}
\title{Sparse canonical covariance analysis
}
\description{
scca is used to perform sparse canonical covariance analysis (SCCA).
scca3 is the extension of scca to address 3 sets of variables on the same set of subjects.
}
\usage{
scca(X,Y,penalty="HL",lamx=c(1,2,3),lamy=c(1,2,3),nc=1,tuning="CV.alt",K=5,seed=NULL,center=TRUE,scale=FALSE)\cr
scca3(X,Y,Z, penalty="HL",lamx=c(1,2,3),lamy=c(1,2,3),lamz=c(1,2,3),nc=1,tuning="CV.alt",K=5,seed=NULL,center=TRUE,scale=FALSE)
}
\arguments{
\tabular{ll}{
X: \tab n-by-p data matrix, where n is the number of subjects and p is the number of variables\cr
Y: \tab n-by-q data matrix, where q is the number of variables\cr
Z: \tab n-by-r data matrix, where r is the number of variables\cr
penalty: \tab "HL" is the unbounded penalty proposed by Lee and Oh (2009). "LASSO" (Tibshirani, 1996), "SCAD" (Fan and Li, 2001)
and "SOFT" (soft thresholding) are also available as other penalty options. Default is "HL". \cr
lamx: \tab A vector specifying grid points of the tuning parameter for X. Default is (1,2,3). \cr
lamy: \tab A vector specifying grid points of the tuning parameter for Y. Default is (1,2,3). \cr
lamz: \tab A vector specifying grid points of the tuning parameter for Z. Default is (1,2,3). \cr
nc: \tab Number of components (canonical vectors). Default is 1. \cr
Tuning: \tab How to find optimal tuning parameters for the sparsity.  If tuning="CV.full", then the tuning parameters are selected
automatically via K-fold cross-validation by using 2-dim'l grid search. If "CV.alt", then a sequential 1-dim'l search method
is applied instead of the 2-dim'l grid search. For scca3, only "CV.alt" is available. Default is "CV.alt". \cr
K: \tab Perform K-fold cross-valiadation. Default is 5. \cr
seed: \tab Seed number for initialization. A random initial point is generated for tuning="CV.alt". \cr
center: \tab The columns of the data matrix are centered to have mean zero. Default is TRUE. \cr
scale: \tab The columns of the data matrix are scaled to have variance 1. Default is FALSE. \cr
}
}
\details{
Sparse CCA uses a random-effect model approach to obtain sparse regression.
This model gives unbounded gains for zero loadings at the origin. Various penalty functions can be adapted as well. 
}
\value{
\tabular{ll}{
A: \tab p-by-nc matrix, k-th colum of A corresponds to k-th pattern (canonical vector) for X \cr
B: \tab q-by-nc matrix, k-th colum of B corresponds to k-th pattern (canonical vector) for Y \cr
U: \tab n-by-nc matrix. k-th column of U corresponds to k-th score associated with k-th pattern for X \cr
V: \tab n-by-nc matrix. k-th column of V corresponds to k-th score associated with k-th pattern for Y \cr
lambda: \tab nc-by-2 matrix. k-th row of lambda corresponds to the optimal tuning parameters for k-th pattern pairs \cr
CR: \tab average cross-validated sample covariance\cr
}
}
\author{
Woojoo Lee, Donghwan Lee, Youngjo Lee and Yudi Pawitan
}
\references{
Lee, W., Lee, D., Lee, Y. and Pawitan, Y. (2011) Sparse Canonical Covariance Analysis for High-throughput Data
}
\examples{

## Example 1
## A very simple simulation example

n<-10; p<-50; q<-20
X = matrix(rnorm(n*p),ncol=p)
Y = matrix(rnorm(n*q),ncol=q)

scca(X,Y)

## Example 2
## Simulation setting in Lee et al.(2011)

n<-50; p<-150; q<-100; r<-120

pcc<-5

sig.e<-0.1
sig.mu<-2 

tx<-rnorm(pcc,1,0.1)
tx.norm<-sqrt(sum(tx^2))
ax<-tx/tx.norm; 
ax<-c(ax,rep(0,p-pcc))

ty<-rnorm(pcc,1,0.1)
ty.norm<-sqrt(sum(ty^2))
by<-ty/ty.norm; 
by<-c(by,rep(0,q-pcc))

tz<-rnorm(pcc,1,0.1)
tz.norm<-sqrt(sum(tz^2))
cz<-tz/tz.norm; 
cz<-c(cz,rep(0,r-pcc))

mu<-rnorm(n,0,sig.mu)

X<-matrix(0,n,p)
Y<-matrix(0,n,q)
Z<-matrix(0,n,r)
for (i in 1:n){ 
    for (j in 1:pcc){
        X[i,j]<-ax[j]*mu[i]+rnorm(1,0,sig.e)
        Y[i,j]<-by[j]*mu[i]+rnorm(1,0,sig.e)
        Z[i,j]<-cz[j]*mu[i]+rnorm(1,0,sig.e)
    }
}
for (i in 1:n){
    for (j in (pcc+1):p){
        X[i,j]<-rnorm(1,0,sig.e)      
    }
}
for (i in 1:n){
    for (j in (pcc+1):q){
        Y[i,j]<-rnorm(1,0,sig.e)      
    }
}

for (i in 1:n){
    for (j in (pcc+1):r){
        Z[i,j]<-rnorm(1,0,sig.e)      
    }
}

scca(X,Y)

scca3(X,Y,Z)

scca(X, Y, penalty="SOFT", lamx=c(0.01,0.05),lamy=c(0.01,0.05), center=FALSE)


}


