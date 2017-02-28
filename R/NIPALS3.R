NIPALS3 <-
function(X,Y,Z){
    eps.cov<-1
    v1<-Y[,1]
    w1<-Z[,1]
    
    
    while (eps.cov>10^-6){
    a1<-t(X)%*%(v1+w1)/as.numeric(t(v1+w1)%*%(v1+w1))
    a1<-a1/sqrt(as.numeric(t(a1)%*%a1))
    u1<-X%*%a1
    
    a1.old<-a1
    
    b1<-t(Y)%*%(w1+u1)/as.numeric(t(w1+u1)%*%(w1+u1))
    b1<-b1/sqrt(as.numeric(t(b1)%*%b1))
    v1<-Y%*%b1
   
    
    c1<-t(Z)%*%(u1+v1)/as.numeric(t(u1+v1)%*%(u1+v1))
    c1<-c1/sqrt(as.numeric(t(c1)%*%c1))
    w1<-Z%*%c1
    
    a1<-t(X)%*%(v1+w1)/as.numeric(t(v1+w1)%*%(v1+w1))
    a1<-a1/sqrt(as.numeric(t(a1)%*%a1))
    
    rho12<-t(u1)%*%v1
    rho23<-t(v1)%*%w1
    rho31<-t(w1)%*%u1
    
    eps.cov<-max(abs(a1-a1.old))
    }
    rho=rho12+rho23+rho31
    return(list(rho=rho,rho12=rho12,rho23=rho23,rho31=rho31,a1=round(a1,digit=4),b1=round(b1,digit=4),c1=round(c1,digit=4)))   
}

