from pylab import *
from HTMLgen import *
from Tables import *
#valArray is a 1d array of scalar values
#fileBase is the baseName for the .png, .ps, and .dat files.
#returns a table containing the plot and links to the .ps and .dat
#produce scalar trace image and ps with data as a function of index
def ProduceTracePicture(valArray,fileBase,hlabel,vlabel,myTitle=''):
    clf()
#    x=fromfunction(lambda i:i,(len(valArray)))

    x=zeros(len(valArray))+0.0
    for i in range(0,len(x)):
        x[i]=i
    plot(x,valArray)
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
    n = len(valArray)
    for i in range(0,len(valArray)):
         asciiFile.write('%20.16e\n' % valArray[i])
    asciiFile.close()

#build table with image, ps, and ascii file in it
    fileTable=NewTable()
    fileTable.width='100%'
    fileTable.body= [[Href(fileBase+".ps",'Postscript'),Href(asciiFileName,'ASCII data')]]
    myTable=NewTable()
    myTable.body=[[myImg]]
    myTable.body.append([fileTable])
    return myTable


#list (over processors) of 1d arrays of scalar data 
def BuildScalarTracePage(valArrayList,baseName,varName,StartCut):
     doc=SimpleDocument()
     row = []
     for i in range(0,len(valArrayList)):
          d = valArrayList[i]
          row.append(ProduceTracePicture(d[StartCut:-1], baseName+"_"+repr(i), 'Blocks', varName,"Proc "+repr(i)))
     picTable = NewTable()
     picTable.width = 350*len(valArrayList)
     picTable.body = [[]]
     picTable.body.append(row)
     doc.append(picTable)
     doc.write(baseName+'.html')
     return baseName+'.html'                    
