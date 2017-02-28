NIPALS.soft3 <-
function(X,Y,Z,lamx,lamy,lamz){
    eps.cov<-1
    n<-dim(X)[1]
    p<-dim(X)[2]
    q<-dim(Y)[2]
    r<-dim(Z)[2]
    v1<-Y[,1]
    w1<-Z[,1]
    res.NIPALS<-NIPALS3(X,Y,Z)
    a1<-res.NIPALS$a1
    b1<-res.NIPALS$b1
    c1<-res.NIPALS$c1
    niter<-0
    
    while ((eps.cov>10^-3)&(niter<100)){
    
    a1.old<-a1
    
    v1w1.norm<-sum((v1+w1)^2)
    Xtv1w1<-crossprod(X, (v1+w1))
    
    if (v1w1.norm>0){ 
    for (i in 1:p){
    a1[i]<-soft((Xtv1w1[i]/v1w1.norm),lamx)
    }
    }
    if (v1w1.norm==0){a1<-rep(0,p)}

    if (sum(a1==0)<p) {a1<-a1/sqrt(sum(a1^2))}
    if (sum(a1==0)==p) {u1<-rep(0,n)}
    u1<-X%*%a1
    
    b1.old<-b1
    
    w1u1.norm<-sum((w1+u1)^2)
    Ytw1u1<-crossprod(Y, (w1+u1))
    if (w1u1.norm>0){
    for (i in 1:q){
    b1[i]<-soft(Ytw1u1[i]/w1u1.norm,lamy)
    }
    }
    if (w1u1.norm==0){b1<-rep(0,q)}

    if (sum(b1==0)<q) {b1<-b1/sqrt(sum(b1^2))}
    if (sum(b1==0)==q) {v1<-rep(0,n)}
    
    v1.old<-v1
    v1<-Y%*%b1


    c1.old<-c1
    
    u1v1.norm<-sum((u1+v1)^2)
    Ztu1v1<-crossprod(Z, (u1+v1)) 
    if (u1v1.norm>0){
    for (i in 1:r){
    c1[i]<-soft(Ztu1v1[i]/u1v1.norm,lamz)
    }
    }
    if (u1v1.norm==0){c1<-rep(0,r)}

    if (sum(c1==0)<r) {c1<-c1/sqrt(sum(c1^2))}
    if (sum(c1==0)==r) {w1<-rep(0,n)}
    
    w1.old<-w1
    w1<-Z%*%c1


    eps.cov<-max(abs(w1-w1.old),abs(v1-v1.old))
    niter<-niter+1
    }
  
    rho12<-t(u1)%*%v1
    rho23<-t(v1)%*%w1
    rho31<-t(w1)%*%u1
    rho<-rho12+rho23+rho31
    return(list(rho=rho,rho12=rho12,rho23=rho23,rho31=rho31,u1=u1,v1=v1,w1=w1,a1=round(a1,digit=4),b1=round(b1,digit=4),c1=round(c1,digit=4),niter=niter)) 
}

