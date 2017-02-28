scca <-
function(X, Y, penalty="HL", lamx=c(1,2,3),lamy=c(1,2,3), nc=1, tuning="CV.alt",K=5, seed=NULL, center=TRUE, scale=FALSE){

    if (penalty=="SOFT" && min(lamx)>=1) {stop("Range of lamx for SOFT should be (0,1)")}
    if (penalty=="SOFT" && min(lamy)>=1) {stop("Range of lamy for SOFT should be (0,1)")}
 
    X<-scale(X,center=center,scale=scale);Y<-scale(Y,center=center,scale=scale)
 
    if (is.null(seed)) {seed<-Sys.time()}             
    if (nrow(X)!=nrow(Y)) {stop("X and Y should have same number of rows")}     

    svd.X<-svd(X) 
    
    n<-dim(X)[1];p<-dim(X)[2]; q<-dim(Y)[2]
    
    U<-matrix(0,n,nc)    
    V<-matrix(0,n,nc)    
    A<-matrix(0,p,nc)    
    B<-matrix(0,q,nc)    
    CR<-rep(0,nc)

    L<-matrix(0,nc,2)        
    
    X.new<-X; Y.new<-Y
    
    for (cnt in 1:nc){
    
        if (tuning[1]=="CV.full") {cv.full<-opt.cv.full(X.new,Y.new,lamx,lamy,K=K,penalty=penalty,seed=seed); olamx<-cv.full$optx; olamy<-cv.full$opty; CR1<-cv.full$avecvcov}
        if (tuning[1]=="CV.alt") {cv.alt<-opt.cv.alt(X.new,Y.new,lamx,lamy,K=K,penalty=penalty,seed=seed); olamx<-cv.alt$optx; olamy<-cv.alt$opty; CR1<-cv.alt$avecvcov}
        
        if (is.numeric(tuning)) {olamx<-tuning[cnt,1];olamy<-tuning[cnt,2]}
        
        
        if (penalty!="SOFT") {sscca.rst<-NIPALS.sparse(X.new,Y.new,olamx, olamy, penalty=penalty)}
        if (penalty=="SOFT") {sscca.rst<-NIPALS.soft(X.new,Y.new,olamx, olamy)}

              
X.new<-X.new- sscca.rst$u1 %*%t(sscca.rst$a1)
        Y.new<-Y.new- sscca.rst$v1 %*%t(sscca.rst$b1)
        

        ytx.new<-t(Y.new)%*%X.new
        L[cnt,1]<-olamx; L[cnt,2]<-olamy

        cat("Computng Component number ",cnt,"\n")
    
        U[,cnt]<-sscca.rst$u1; V[,cnt]<-sscca.rst$v1
        A[,cnt]<-sscca.rst$a1; B[,cnt]<-sscca.rst$b1
CR[cnt]<-CR1

    }

    
return(list(U=U, V=V, A=A, B=B, lambda=L,CR=CR))
}

