opt.cv.alt <-
function(X,Y,K, lamx, lamy, penalty,seed){  

    row.x<-dim(X)[1]
    seq.x<-lamx
    seq.y<-lamy

    i.x<-length(seq.x)
    i.y<-length(seq.y)
    cvcov<-rep(0,K)
    cvcovdiff<-rep(0,K)
    avecvcov<-matrix(0,i.x,i.y)
    avecvcovdiff<-matrix(0,i.x,i.y)
       
    set.seed(seed)

    cv.set<-split(sample(1:row.x),rep(1:K,length=row.x))
    
    opt.i<-sample(c(1:i.x),1) 
    
    opt.conv<-1;max.old<-0;opt.iter<-0
    
    while(opt.conv>1e-02 & opt.iter<20){
    
        for (jj in 1:i.y){
    
            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[opt.i],seq.y[jj], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[opt.i],seq.y[jj])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            
            cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
            cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
            }
            
            avecvcov[opt.i,jj]<-sum(abs(cvcov))/K
            avecvcovdiff[opt.i,jj]<-sum(abs(cvcovdiff))/K
        }
        
    opt.j<- which.max(avecvcov[opt.i,])
    
            for (ii in 1:i.x){
    
            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse(XX,YY,seq.x[ii],seq.y[opt.j], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft(XX,YY,seq.x[ii],seq.y[opt.j])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            
            cvcov[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))
            cvcovdiff[kk]<-abs(t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))-abs(resXY$rho1)
            }
            
            avecvcov[ii,opt.j]<-sum(abs(cvcov))/K
            avecvcovdiff[ii,opt.j]<-sum(abs(cvcovdiff))/K
        }
    opt.i<- which.max(avecvcov[,opt.j])    
    
    max.new<-max(avecvcov[opt.i,opt.j])
    opt.conv<-abs(max.new-max.old)  
    max.old<-max.new  
    opt.iter<-opt.iter+1
    }
        

    r.index<-opt.i
    c.index<-opt.j
    
    return(list(optx=seq.x[r.index],opty=seq.y[c.index],avecvcov=avecvcov[r.index,c.index]))
}

