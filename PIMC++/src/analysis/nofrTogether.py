#import numarray
import numpy
import pickle
import math
from pylab import *


def ProduceCorrelationPicture(x,y,fits,fileBase,hlabel,vlabel,cutoff,both):
     clf()
     #if (IsMonotonic(x)):
     #     plot(x, y)
     #else:
     numLines=len(y)
     xp = numpy.linspace(0, cutoff, 100)
     for line in range(0,numLines):
       if both:
         plot(x[line], y[line], '+')
       fitLine = map(lambda x: fits[line][0]*x + fits[line][1],xp.tolist())
       plot(xp.tolist(), fitLine)
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
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")

def ProduceCorrelationPictureLogLog(x,y,fits,fileBase,hlabel,vlabel,cutoff):
     clf()
     #if (IsMonotonic(x)):
     #     plot(x, y)
     #else:
     numLines=len(y)
     xp = numpy.linspace(0.001, cutoff, 100)
     for line in range(0,numLines):
       loglog(x[line], y[line], '+', basex=10)
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
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")

##Produce Image
f = open('nofrlineDump', 'r')
nLines = int(sys.argv[1])
x=[]
y=[]
for i in range(0,nLines):
  x.append(pickle.load(f))
  y.append(pickle.load(f))
f.close()
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
myImgTable=ProduceCorrelationPicture(x,y,fits,baseName,hlabel,vlabel,box,1)


fits=[]
xCut=[]
yCut=[]
logx=[]
logy=[]
print "Making Figure 2"
for i in range(0,nLines):
  xCut.append(log(x[i][70:100])/log(10))
  yCut.append(log(y[i][70:100])/log(10))
  fits.append(numpy.polyfit(xCut[i],yCut[i],1))
  logx.append(log(x[i])/log(10))
  logy.append(log(y[i])/log(10))
myImgTable2=ProduceCorrelationPicture(logx,logy,fits,baseName+"loglog","log(r)","log(n(r))",log(box)/log(10),0)

print "Making Figure 3"
myImgTable3=ProduceCorrelationPicture(logx,logy,fits,baseName+"loglog2","log(r)","log(n(r))",log(box)/log(10),1)

print "Making Figure 4"
myImgTable4=ProduceCorrelationPictureLogLog(x,y,fits,baseName+"loglog3",hlabel,vlabel,box)
