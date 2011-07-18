import numpy
import stats
import math
from pylab import *

def ProcessData(InputFiles, Observable):

  numLines = len(InputFiles)
  FileList = []
  Base = []
  RawData = []
  Data = []
  x = []
  y = []
  err = []
  key = []
  for i in range(0, numLines):
    FileList.append(open(InputFiles[i]).readlines())

    j = 0
    for line in FileList[i]:
      FileList[i][j] = line.strip()
      j += 1
    print FileList[i]

    key.append(str(FileList[i][0]))
    FileList[i].remove(FileList[i][0])

    Base.append([])
    RawData.append([])
    j = 0
    for File in FileList[i]:
      Base[i].append(File.replace('.0.h5',''))
      RawData[i].append(open(Base[i][j] + '/' + Observable + '.dat').readlines())
      j += 1
    # print RawData[i]

    Data.append([])
    for entry in RawData[i]:
      Data[i].append(entry[0].split())
    # print Data[i]

    x.append([])
    y.append([])
    err.append([])
    for j in range(0,len(Data[i])):
      x[i].append(float(Data[i][j][0]))
      y[i].append(float(Data[i][j][1]))
      err[i].append(float(Data[i][j][2]))

    print x[i], y[i], err[i], key[i]

  return x, y, err, key

def PlotData(x,y,err,key,fileBase,hlabel,vlabel):
  clf()
  numLines=len(y)
  for line in range(0,numLines):
    errorbar(x[line], y[line], xerr=None, yerr=err[line], label=key[line])
  h1=xlabel(hlabel)
  setp(h1,"FontSize",20)
  v1=ylabel(vlabel)
  setp(v1,"FontSize",20)
  labels = get(gca(), 'xticklabels')
  setp(labels, 'fontsize', 16)
  labels = get(gca(), 'yticklabels')
  setp(labels, 'fontsize', 16)
  currAxis=list(axis())
  axis(currAxis)
  grid(linestyle=':', linewidth=1)
  legend(markerscale=2)
  savefig(fileBase+".png",dpi=60)
  savefig(fileBase+".ps")

Observable = sys.argv[1]

InputFiles = []
for i in range(0, len(sys.argv)-2):
  InputFiles.append(sys.argv[i+2])
# print InputFiles

(x,y,err,key) = ProcessData(InputFiles, Observable)
PlotData(x,y,err,key,Observable,'Ry',Observable)

