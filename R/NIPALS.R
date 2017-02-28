NIPALS <-
function(X,Y){
    eps.cov<-1
    v1<-Y[,1]

    while (eps.cov>10^-6){
        a1<-crossprod(X,v1)/sum(v1^2)
        a1<-a1/sqrt(sum(a1^2))
        u1<-X%*%a1
        
        b1<-crossprod(Y,u1)/sum(u1^2)
        b1<-b1/sqrt(sum(b1^2))
        v1.old<-v1
        v1<-Y%*%b1
        eps.cov<-max(abs(v1-v1.old))
    }
    rho1<-crossprod(u1, v1)
    return(list(rho1=rho1,u1=u1,v1=v1,a1=round(a1,digit=4),b1=round(b1,digit=4)))   
}

