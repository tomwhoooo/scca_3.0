scca3 <-
function(X, Y, Z, penalty="HL", lamx=c(1,2,3), lamy=c(1,2,3),lamz=c(1,2,3), nc=1, tuning="CV.alt",K=5, seed=NULL,center=TRUE,scale=FALSE){


    if (penalty=="SOFT" && min(lamx)>=1) {stop("Range of lamx for SOFT should be (0,1)")}
    if (penalty=="SOFT" && min(lamy)>=1) {stop("Range of lamy for SOFT should be (0,1)")}
    if (penalty=="SOFT" && min(lamz)>=1) {stop("Range of lamz for SOFT should be (0,1)")}
   
    X<-scale(X,center=center,scale=scale);Y<-scale(Y,center=center,scale=scale);Z<-scale(Z,center=center,scale=scale)
 
    if (is.null(seed)) {seed<-Sys.time()}  
            
    if (nrow(X)!=nrow(Y)) {stop("X and Y should have same number of rows")}     
    if (nrow(Y)!=nrow(Z)) {stop("Y and Z should have same number of rows")} 
    if (nrow(X)!=nrow(Z)) {stop("X and Z should have same number of rows")}   
  
    svd.X<-svd(X)
       
    if (is.null(nc)) {nc<-1}
    n<-dim(X)[1]; p<-dim(X)[2]; q<-dim(Y)[2]; r<-dim(Z)[2]

    U<-matrix(0,n,nc)    
    V<-matrix(0,n,nc)    
    W<-matrix(0,n,nc)    
    A<-matrix(0,p,nc)    
    B<-matrix(0,q,nc)    
    C<-matrix(0,r,nc)
    CR<-rep(0,nc)

    L<-matrix(0,nc,3)         
    
    X.new<-X; Y.new<-Y; Z.new<-Z

    for (cnt in 1:nc){
    
        if (tuning[1]=="CV.alt") {cv.alt<-opt.cv.alt3(X.new,Y.new,Z.new,lamx,lamy,lamz,K=5,penalty=penalty,seed=seed); olamx<-cv.alt$optx; olamy<-cv.alt$opty; olamz<-cv.alt$optz; CR1<-cv.alt$avecvcov}
        
        if (is.numeric(tuning)) {olamx<-tuning[cnt,1];olamy<-tuning[cnt,2];olamz<-tuning[cnt,3]}
        
        
        if (penalty!="SOFT") {sscca.rst<-NIPALS.sparse3(X.new,Y.new,Z.new, olamx, olamy, olamz, penalty=penalty)}
        if (penalty=="SOFT") {sscca.rst<-NIPALS.soft3(X.new,Y.new,Z.new, olamx, olamy, olamz)}

        
X.new<-X.new- sscca.rst$u1 %*%t(sscca.rst$a1)
        Y.new<-Y.new- sscca.rst$v1 %*%t(sscca.rst$b1)
        Z.new<-Z.new- sscca.rst$w1 %*%t(sscca.rst$c1)

        ytx.new<-t(Y.new)%*%X.new

        L[cnt,1]<-olamx; L[cnt,2]<-olamy; L[cnt,3]<-olamz

        cat("Computng Component number ",cnt,"\n")
    
        U[,cnt]<-sscca.rst$u1; V[,cnt]<-sscca.rst$v1; W[,cnt]<-sscca.rst$w1
        A[,cnt]<-sscca.rst$a1; B[,cnt]<-sscca.rst$b1; C[,cnt]<-sscca.rst$c1
CR[cnt]<-CR1

    }

    
return(list(U=U, V=V, W=W,A=A, B=B, C=C,lambda=L,CR=CR))
}

