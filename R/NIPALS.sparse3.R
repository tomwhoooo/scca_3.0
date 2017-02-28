NIPALS.sparse3 <-
function(X,Y,Z,lamx,lamy,lamz, penalty){
    eps.cov<-1
    p<-dim(X)[2]
    q<-dim(Y)[2]
    r<-dim(Z)[2]
    v1<-Y[,1]
    w1<-Z[,1]

    res0<-NIPALS3(X,Y,Z)
    a1<-res0$a1
    b1<-res0$b1
    c1<-res0$c1
    niter<-0
    a<-3.7
    w<-30
    if (lamx==0) {sigx<-median(sd(X))}
    if (lamy==0) {sigy<-median(sd(Y))}
    if (lamz==0) {sigz<-median(sd(Z))}

    if (lamx!=0) {sigx<-sd(res0$a1)/sqrt(2)}
    if (lamy!=0) {sigy<-sd(res0$b1)/sqrt(2)}
    if (lamz!=0) {sigz<-sd(res0$c1)/sqrt(2)} 
    
    if (lamx==0 & penalty=="SCAD") {lamx<-1e-10}
    if (lamy==0 & penalty=="SCAD") {lamy<-1e-10}
    if (lamz==0 & penalty=="SCAD") {lamz<-1e-10}
    
    while ((eps.cov>1e-03)&(niter<100)){
    eps.cov.a1<-1
    nitera<-0
    while ((eps.cov.a1>1e-03)&(nitera<50)){
    
        if (penalty=="LASSO") {uux<-abs(a1)+1e-08}
        if (penalty=="SCAD") {uux<-(abs(a1)+1e-08)/ (as.numeric(abs(a1)<=lamx)+as.numeric(abs(a1)>lamx)*(a*lamx-abs(a1))*as.numeric(a*lamx>abs(a1))/((a-1)*lamx))}
        if (penalty=="HL"){    kax <- sqrt(4*a1^2/(w*sigx^2)+((2/w)-1)^2) ;    uux <- 0.25*w*(((2/w)-1)+kax)+1e-08 }
        
        WWx<-as.vector(1/uux)
        a1.old<-a1
        
        v1w1.norm<-sum((v1+w1)^2)
        Xtv1w1<-crossprod(X, (v1+w1))
        
        for (i in 1:p){
            a1[i]<-(1/(v1w1.norm+lamx*WWx[i]))*Xtv1w1[i]
        }
        eps.cov.a1<-max(abs(a1-a1.old))
    }
    a1.norm<-sum(a1^2)
    a1<-a1/sqrt(a1.norm) 

    if (nitera>1) {u1.old<-u1}
    u1<-X%*%a1
    
    
    eps.cov.b1<-1
    niterb<-0
    while ((eps.cov.b1>1e-03)&(niterb<50) ){
    
        if (penalty=="LASSO") {uuy<-abs(b1)+1e-08}
        if (penalty=="SCAD") {uuy<-(abs(b1)+1e-08)/ (as.numeric(abs(b1)<=lamy)+as.numeric(abs(b1)>lamy)*(a*lamy-abs(b1))*as.numeric(a*lamy>abs(b1))/((a-1)*lamy))}
        if (penalty=="HL"){    kay <- sqrt(4*b1^2/(w*sigy^2)+((2/w)-1)^2) ;    uuy <- 0.25*w*(((2/w)-1)+kay)+1e-08 }
        
        WWy<-as.vector(1/uuy)
        b1.old<-b1
        
        w1u1.norm<-sum((w1+u1)^2)
        Ytw1u1<-crossprod(Y, (w1+u1))
            for (i in 1:q){
                b1[i]<-(1/(w1u1.norm+lamy*WWy[i]))*Ytw1u1[i]
            }
        eps.cov.b1<-max(abs(b1-b1.old))
        }
        

    b1.norm<-sum(b1^2)
    b1<-b1/sqrt(b1.norm) 
    
    v1.old<-v1
    v1<-Y%*%b1


    eps.cov.c1<-1
    niterc<-0
    while ((eps.cov.c1>1e-03)&(niterc<50) ){
    
        if (penalty=="LASSO") {uuz<-abs(c1)+1e-08}
        if (penalty=="SCAD") {uuz<-(abs(c1)+1e-08)/ (as.numeric(abs(c1)<=lamz)+as.numeric(abs(c1)>lamz)*(a*lamz-abs(c1))*as.numeric(a*lamz>abs(c1))/((a-1)*lamz))}
        if (penalty=="HL"){ kaz <- sqrt(4*c1^2/(w*sigz^2)+((2/w)-1)^2) ;    uuz <- 0.25*w*(((2/w)-1)+kaz)+1e-08 }
        
        WWz<-as.vector(1/uuz)
        c1.old<-c1
        
        u1v1.norm<-sum((u1+v1)^2)
        Ztu1v1<-crossprod(Z, (u1+v1))
            for (i in 1:r){
                c1[i]<-(1/(u1v1.norm+lamz*WWz[i]))*Ztu1v1[i]
            }
        eps.cov.c1<-max(abs(c1-c1.old))
        }
        

    c1.norm<-sum(c1^2)
    c1<-c1/sqrt(c1.norm) 
    
    w1.old<-w1
    w1<-Z%*%c1



    eps.cov<-max(abs(w1-w1.old))
    niter<-niter+1
    if ((sum(round(a1,digit=4)==0)>=(p-2))&(sum(round(b1,digit=4)==0)>=(q-2))) {break}

    }
        
    a1<-a1/sqrt(sum(a1^2))
    b1<-b1/sqrt(sum(b1^2))
    c1<-c1/sqrt(sum(c1^2))   
    u1<-X%*%a1  
    v1<-Y%*%b1
    w1<-Z%*%c1
    rho12<-crossprod(u1, v1)
    rho23<-crossprod(v1, w1)
    rho31<-crossprod(w1, u1)
    rho<-rho12+rho23+rho31
    if (niter>=100) {print("No Convergence")}
    return(list(rho=rho,rho12=rho12,rho23=rho23,rho31=rho31,u1=u1,v1=v1,w1=w1,a1=round(a1,digit=4),b1=round(b1,digit=4),c1=round(c1,digit=4),WWx=WWx,WWy=WWy,WWz=WWz,niter=niter))   
    
}

