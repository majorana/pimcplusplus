from numarray import *
def c(i,x,mean,var):
    N=len(x)
    corr=1.0/var*1.0/(N-i)*sum((x[0:N-i]-mean)*(x[i:N]-mean))
    return corr
                         
def Stats(x):
    N=len(x)
    mean=sum(x)/(N+0.0)
    xSquared=x*x
    var=sum(xSquared)/(N+0.0)-mean*mean
    i=0          
    tempC=0.5
    kappa=0.0
    while (tempC>0 and i<N):
        kappa=kappa+2.0*tempC
        i=i+1
        tempC=c(i,x,mean,var)
    Neff=(N+0.0)/(kappa+0.0)
    error=sqrt(var/Neff)
    return (mean,var,error,kappa)

