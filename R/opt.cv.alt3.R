opt.cv.alt3 <-
function(X,Y,Z,K,lamx,lamy,lamz,penalty,seed){  

    row.x<-dim(X)[1]
    seq.x<-lamx; seq.y<-lamy; seq.z<-lamz
 
    i.x<-length(seq.x)
    i.y<-length(seq.y)
    i.z<-length(seq.z) 
    cvcov<-rep(0,K)
    cvcovdiff<-rep(0,K)
    avecvcov<-array(0,c(i.x,i.y,i.z))
    avecvcovdiff<-array(0,c(i.x,i.y,i.z))
    
    set.seed(seed)    

    cv.set<-split(sample(1:row.x),rep(1:K,length=row.x)) 
        
    opt.i<-sample(c(1:i.x),1) 
    opt.k<-sample(c(1:i.z),1) 

    opt.conv<-1;max.old<-0;opt.iter<-0
    
    while(opt.conv>1e-02 & opt.iter<20){
    
        for (jj in 1:i.y){

            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            ZZ<-Z[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[jj],seq.z[opt.k], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[jj],seq.z[opt.k], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[jj],seq.z[opt.k], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft3(XX,YY,ZZ,seq.x[opt.i],seq.y[jj],seq.z[opt.k])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            rc1<-resXY$c1

            cvcov[kk]<-abs((t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))+(t(Y[cv.set[[kk]],]%*%rb1)%*%(Z[cv.set[[kk]],]%*%rc1))+(t(Z[cv.set[[kk]],]%*%rc1)%*%(X[cv.set[[kk]],]%*%ra1)))
            }
            
            avecvcov[opt.i,jj,opt.k]<-sum(abs(cvcov))/K
          
       }
        
    opt.j<- which.max(avecvcov[opt.i,,opt.k])
    
            for (ii in 1:i.x){
    
            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            ZZ<-Z[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[ii],seq.y[opt.j],seq.z[opt.k], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[ii],seq.y[opt.j],seq.z[opt.k], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[ii],seq.y[opt.j],seq.z[opt.k], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft3(XX,YY,ZZ,seq.x[ii],seq.y[opt.j],seq.z[opt.k])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            rc1<-resXY$c1

            cvcov[kk]<-abs((t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))+(t(Y[cv.set[[kk]],]%*%rb1)%*%(Z[cv.set[[kk]],]%*%rc1))+(t(Z[cv.set[[kk]],]%*%rc1)%*%(X[cv.set[[kk]],]%*%ra1)))
            }
            
            avecvcov[ii,opt.j,opt.k]<-sum(abs(cvcov))/K
        }
    opt.i<- which.max(avecvcov[,opt.j,opt.k])    
    
    
   for (mm in 1:i.z){
    
            for (kk in 1:K){
            XX<-X[-cv.set[[kk]],]
            YY<-Y[-cv.set[[kk]],]
            ZZ<-Z[-cv.set[[kk]],]
            if (penalty=="LASSO") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[opt.j],seq.z[mm], "LASSO")}
            if (penalty=="SCAD") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[opt.j],seq.z[mm], "SCAD")}
            if (penalty=="HL") {resXY<-NIPALS.sparse3(XX,YY,ZZ,seq.x[opt.i],seq.y[opt.j],seq.z[mm], "HL")}
            if (penalty=="SOFT") {resXY<-NIPALS.soft3(XX,YY,ZZ,seq.x[opt.i],seq.y[opt.j],seq.z[mm])}
            
            ra1<-resXY$a1
            rb1<-resXY$b1
            rc1<-resXY$c1

            cvcov[kk]<-abs((t(X[cv.set[[kk]],]%*%ra1)%*%(Y[cv.set[[kk]],]%*%rb1))+(t(Y[cv.set[[kk]],]%*%rb1)%*%(Z[cv.set[[kk]],]%*%rc1))+(t(Z[cv.set[[kk]],]%*%rc1)%*%(X[cv.set[[kk]],]%*%ra1)))
            }
            
            avecvcov[opt.i,opt.j,mm]<-sum(abs(cvcov))/K
        }

   opt.k<- which.max(avecvcov[opt.i,opt.j,]) 




    max.new<-max(avecvcov[opt.i,opt.j,opt.k])
    opt.conv<-abs(max.new-max.old)  
    max.old<-max.new  
    opt.iter<-opt.iter+1
    }
        

    r.index<-opt.i
    c.index<-opt.j
    h.index<-opt.k
    
    return(list(optx=seq.x[r.index],opty=seq.y[c.index],optz=seq.z[h.index],avecvcov=avecvcov[r.index,c.index,h.index]))
}

