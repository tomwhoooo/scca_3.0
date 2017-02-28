NIPALS.soft <-
function(X,Y,lamx,lamy){
    eps.cov<-1
    p<-dim(X)[2]
    q<-dim(Y)[2]
    v1<-Y[,1]
    res.NIPALS<-NIPALS(X,Y)
    a1<-res.NIPALS$a1
    b1<-res.NIPALS$b1
    niter<-0
    
    while ((eps.cov>10^-3)&(niter<100)){
    
    a1.old<-a1
    
    v1.norm<-sum(v1^2)
    Xtv1<-crossprod(X, v1)
    
    for (i in 1:p){
    a1[i]<-soft((Xtv1[i]/v1.norm),lamx)
    }
    
    if (sum(a1==0)<p) {a1<-a1/sqrt(sum(a1^2))}
    if (sum(a1==0)==p) {u1<-rep(0,n);b1<-rep(0,q);v1<-rep(0,n);break}
    u1<-X%*%a1
    
    b1.old<-b1
    
    u1.norm<-sum(u1^2)
    Ytu1<-crossprod(Y, u1)
    for (i in 1:q){
    b1[i]<-soft(Ytu1[i]/u1.norm,lamy)
    }
    
    if (sum(b1==0)<q) {b1<-b1/sqrt(sum(b1^2))}
    if (sum(b1==0)==q) {u1<-rep(0,n);a1<-rep(0,p);v1<-rep(0,n);break}
    
    v1.old<-v1
    v1<-Y%*%b1
    eps.cov<-max(abs(v1-v1.old))
    niter<-niter+1
    }

    rho1<-t(u1)%*%v1
    return(list(rho1=rho1,u1=u1,v1=v1,a1=round(a1,digit=4),b1=round(b1,digit=4),niter=niter))   
}

