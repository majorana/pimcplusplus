#import numarray
import numpy
import pickle
import math
from pylab import *


def ProduceCorrelationPicture(x,y,err,fits,fileBase,hlabel,vlabel,cutoff,both,key):
     clf()
     #if (IsMonotonic(x)):
     #     plot(x, y)
     #else:
     numLines=len(y)
     xp = numpy.linspace(0, cutoff, 100)
     for line in range(0,numLines):
       if both:
         errorbar(x[line], y[line], xerr=None, yerr=err[line], label=key[line])
       fitLine = map(lambda x: fits[line][0]*x + fits[line][1],xp.tolist())
       plot(xp.tolist(), fitLine, label=key[line])
     h1=xlabel(hlabel)
     setp(h1,"FontSize",20)
     v1=ylabel(vlabel)
     setp(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     setp(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     setp(labels, 'fontsize', 16)
     currAxis=list(axis())
     #currAxis[1]=cutoff
     axis(currAxis)
     grid(linestyle=':', linewidth=1)
     #legend(markerscale=2)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")

def ProduceCorrelationPictureLogLog(x,y,fits,fileBase,hlabel,vlabel,cutoff,key):
     clf()
     #if (IsMonotonic(x)):
     #     plot(x, y)
     #else:
     numLines=len(y)
     xp = numpy.linspace(0.001, cutoff, 100)
     for line in range(0,numLines):
       loglog(x[line], y[line], '+', label=key[line], basex=10)
       fitLine = map(lambda x: (10**(fits[line][1]/fits[line][0])) * (x**fits[line][0]),xp.tolist())
       #loglog(xp.tolist(), fitLine, basex=10)
     h1=xlabel(hlabel)
     setp(h1,"FontSize",20)
     v1=ylabel(vlabel)
     setp(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     setp(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     setp(labels, 'fontsize', 16)
     currAxis=list(axis())
     #currAxis[1]=cutoff
     axis(currAxis)
     ylim(0.001,1.001)
     #legend()
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")

##Produce Image
f1 = open('nofrlineDump', 'r')
f2 = open('list.txt','r')
nLines = int(sys.argv[1])
x=[]
y=[]
err=[]
key=[]
for line in f2:
  key.append(line)
for i in range(0,nLines):
  x.append(pickle.load(f1))
  y.append(pickle.load(f1))
  err.append(pickle.load(f1))
f1.close()
baseName="nofrTogether"
hlabel="r"
vlabel="n(r)"
box=30
#for i in range(0,len(x[0])):
#  print x[0][i], i
fits=[]
xCut=[]
yCut=[]
for i in range(0,nLines):
  xCut.append(x[i][70:100])
  yCut.append(y[i][70:100])
  fits.append(numpy.polyfit(xCut[i],yCut[i],1))
print "Making Figure 1"
myImgTable=ProduceCorrelationPicture(x,y,err,fits,baseName,hlabel,vlabel,box,1,key)

logerr=[]
fits=[]
xCut=[]
yCut=[]
logx=[]
logy=[]
print "Making Figure 2"
for i in range(0,nLines):
  y1=asarray(y[i])+asarray(err[i])
  y2=asarray(y[i])-asarray(err[i])

  for j in range(0,len(y1)):
    if y1[j]!=0:
      y1[j]=math.log(y1[j])/math.log(10)
    else:
      y1[j]=0
  for j in range(0,len(y2)):
    if y2[j]>=0:
      y2[j]=math.log(y2[j])/math.log(10)
    else:
      y2[j]=0
  logerr.append(((y1-y2)/2).tolist())

  xCut.append(log(x[i][70:100])/log(10))
  yCut.append(log(y[i][70:100])/log(10))
  fits.append(numpy.polyfit(xCut[i],yCut[i],1))
  logx.append(log(x[i])/log(10))
  logy.append(log(y[i])/log(10))
myImgTable2=ProduceCorrelationPicture(logx,logy,logerr,fits,baseName+"loglog","log(r)","log(n(r))",log(box)/log(10),0,key)

print "Making Figure 3"
myImgTable3=ProduceCorrelationPicture(logx,logy,logerr,fits,baseName+"loglog2","log(r)","log(n(r))",log(box)/log(10),1,key)

print "Making Figure 4"
myImgTable4=ProduceCorrelationPictureLogLog(x,y,fits,baseName+"loglog3",hlabel,vlabel,box,key)
