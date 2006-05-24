import numarray
import IO
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats

def ProcessLangevin (infiles, summaryDoc, detailedDoc, StartCut, beta):
    V = infiles.ReadVar("V")[0]
    bands = infiles.ReadVar("BandEnergies")
    numPtcls = V.shape[1]
    KE = 0.5*numarray.sum(numarray.sum(V*V, 2),1)
    numSteps = len(KE)
    step = arange(0,len(KE))
    expectedKE = 1.5*numPtcls/beta * (step+1)/(step+1)
    baseName  = "LangevinKE"
    imageName = baseName + ".png"
    epsName   = baseName + ".eps"
    
    plot (step, KE, step, expectedKE)
    savefig (imageName, dpi=60)
    savefig (epsName,   dpi=60)
    img = Image(imageName)
    summaryDoc.append(Heading(2, "Langevin dynamics"))
    summaryDoc.append(img)

#    for i in range(0,len(bands)):
#        bands[i] = numarray.transpose(bands[i])
    baseName = "BandEnergies"
    bandPage = BuildScalarTracePage(bands, baseName, "BandEnergies", StartCut)
    myFrame=IFrame()
    myFrame.src=bandPage
    myFrame.width="100%"
    myFrame.height="375"
    detailedDoc.append(myFrame)
    
    
    
