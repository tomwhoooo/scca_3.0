opt.cv.full <-
function(X,Y,K, lamx,lamy,penalty,seed){  

    row.x<-dim(X)[1]
    seq.x<-lamx;  seq.y<-lamy 

    i.x<-length(seq.x)
    i.y<-length(seq.y)
    cvcov<-rep(0,K)
    cvcovdiff<-rep(0,K)
    avecvcov<-matrix(0,i.x,i.y)
    avecvcovdiff<-matrix(0,i.x,i.y)
    
    set.seed(seed)

    cv.set<-split(sample(1:row.x),rep(1:K,length=row.x))

    
    for (ii in 1:i.x){
        for (jj in 1:i.y){
    
            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[jj], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[ii],seq.y[jj])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            
            cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
            cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
            }
            
            avecvcov[ii,jj]<-sum(abs(cvcov))/K
            avecvcovdiff[ii,jj]<-sum(abs(cvcovdiff))/K
        }
    }
    r.index<-(which.max(avecvcov)%%i.x)*((which.max(avecvcov)%%i.x)>0)+i.x*((which.max(avecvcov)%%i.x)==0)
    c.index<-(floor(which.max(avecvcov)/i.x)+1)*((which.max(avecvcov)%%i.x)>0)+(floor(which.max(avecvcov)/i.x))*((which.max(avecvcov)%%i.x)==0)
    
    r.indexw<-(which.min(avecvcovdiff)%%i.x)*((which.min(avecvcovdiff)%%i.x)>0)+i.x*((which.min(avecvcovdiff)%%i.x)==0)
    c.indexw<-(floor(which.min(avecvcovdiff)/i.x)+1)*((which.min(avecvcovdiff)%%i.x)>0)+(floor(which.min(avecvcovdiff)/i.x))*((which.min(avecvcovdiff)%%i.x)==0)
    return(list(optx=seq.x[r.index],opty=seq.y[c.index],optxw=seq.x[r.indexw],optyw=seq.y[c.indexw],avecvcov=avecvcov[r.index,c.index],avecvcovdiff=avecvcovdiff))
}

