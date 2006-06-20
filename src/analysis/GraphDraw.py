from pylab import *
from HTMLgen import *
from Tables import *

def IsMonotonic (x):
     isMono = True
#     for i in range(0,x.size()-2):
     for i in range(0,len(x)-2):
          isMono = isMono and (x[i+1] > x[i])
     return isMono

def ProduceCorrelationPicture(x,y,fileBase,hlabel,vlabel):
     clf()
     #if (IsMonotonic(x)):
     #     plot(x, y)
     #else:
     plot(x, y, 'o')
     h1=xlabel(hlabel)
     setp(h1,"FontSize",20)
     v1=ylabel(vlabel)
     setp(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     setp(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     setp(labels, 'fontsize', 16)
     currAxis=axis()
##     currAxis[0]=cutoff
     axis(currAxis)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")
     myImg=Image(fileBase+".png")
     return myImg



def LongRangeImage(basename,r,long,short,myTitle,labelY):
     clf()
     hold ("off")
     l1 = plot (r, long)
     a = axis()
     hold ("on")
     l2 = plot (r[1:-1], long[1:-1]+short[1:-1], 'r-')
     setp (l1, 'linewidth', 2);
     setp (l2, 'linewidth', 2);
     h1 = xlabel("r")
     axis(a)
     setp (h1, "FontSize", 20)
     h2 = ylabel (labelY)
     setp (h2, "FontSize", 20)
     labels = get(gca(), 'xticklabels')
     setp(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     setp(labels, 'fontsize', 16)
     h3 = legend ('Ulong')
     h4 = title (myTitle)
     setp (h4, "FontSize", 20)
     savefig(basename+".png",dpi=60)
     return Image(basename+".png")


def ProduceTracePicture(data,fileBase,hlabel,vlabel,myTitle=''):
#produce scalar trace image and ps with data as a function of index
    clf()
    x=fromfunction(lambda i:i,(len(data),))
    plot(x,data)
    h1=xlabel(hlabel)
    setp(h1,"FontSize",20)
    v1=ylabel(vlabel)
    setp(v1,"FontSize",20)
    if len(myTitle) != 0:
         t1=title(myTitle)
         setp(t1,"FontSize",20)
    labels = get(gca(), 'xticklabels')
    setp(labels, 'fontsize', 16)
    labels = get(gca(), 'yticklabels')
    setp(labels, 'fontsize', 16)
    savefig(fileBase+".png",dpi=45)
    savefig(fileBase+".ps")
    myImg=Image(fileBase+".png")

#produce asciiFile
    asciiFileName = fileBase + '.dat'
    asciiFile = open (asciiFileName, "w")
    n = len(data)
    for i in range(0,len(data)):
         asciiFile.write('%20.16e\n' % data[i])
    asciiFile.close()

#build table with image, ps, and ascii file in it
    fileTable=NewTable()
    fileTable.width='100%'
    fileTable.body= [[Href(fileBase+".ps",'Postscript'),Href(asciiFileName,'ASCII data')]]
    myTable=NewTable()
    myTable.body=[[myImg]]
    myTable.body.append([fileTable])
    return myTable


