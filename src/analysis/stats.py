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
    tempC=1.0
    kappa=0.0
    while (tempC>0 and i<N):
        kappa=kappa+tempC
        i=i+1
        tempC=c(i,x,mean,var)
    kappa=2.0*kappa+1.0
    Neff=(N+0.0)/(kappa+0.0)
    error=sqrt(var/Neff)
    return (mean,var,error,kappa)

    
        
a=zeros(5)
a[0]=1.
a[1]=3.
a[2]=2.
a[4]=1.
print Stats(a)
