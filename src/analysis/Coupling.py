import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import numarray.fft
import stats
import math



def ProcessCoupling(infiles,summaryDoc,detailedDoc,StartCut):
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("y")
    infoArray=[]
    for proc in range(0,len(data)):
         (mean,var,error,kappa)=stats.Stats(data[proc][StartCut:-1,proc])
         (mean2,var2,error2,kappa2)=stats.Stats(data[proc][StartCut:-1,proc+1])
         infoArray.append((mean,mean2,var,var2,error,error2,kappa,kappa2))
    print mean,var,error
    print mean2,var2,error2
    print mean+mean2
    totalProd=1.0
    totalError2=0.0
    for counter in range(20,len(infoArray)):
         (mean,mean2,_,_,e1,e2,_,_)=infoArray[counter]
         if (not(mean==0 or mean2==0)):
             totalProd=totalProd*(mean/mean2)
         if (not(mean==0 or mean2==0)):
             totalError2=totalError2+(e1/mean)**2+(e2/mean2)**2
         if (mean==0 or mean2==0):
             print "A mean is 0"
         else:
            print mean/mean2, math.sqrt((e1/mean)**2+(e2/mean2)**2),totalProd,math.sqrt(totalError2)*totalProd

    print totalProd
    print math.sqrt(totalError2)*totalProd
    return 0
