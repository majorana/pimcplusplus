#!/usr/bin/python
import tables
import sys
import stats
#import Numeric
h5file = tables.openFile(sys.argv[1])
#data = h5file.root.Paths_2.Path.read()
#names=h5file.root.Paths_2.SpeciesNames.read()

E = h5file.root.Energies_0.TotalEnergy.read()
PE= h5file.root.Energies_0.PotentialEnergy.read()
SE= h5file.root.Energies_0.SpringEnergy.read()
UE=h5file.root.Energies_0.DBetaEnergy.read()
(mean,var,error,kappa)= stats.Stats(E)
(meanP,varP,errorP,kappaP)=stats.Stats(PE)
(meanS,varS,errorS,kappaS)=stats.Stats(SE)
(meanU,varU,errorU,kappaU)=stats.Stats(UE)
print 'Mean:  ', mean
print 'Var:   ', var
print 'Error: ', error
print 'kappa  ', kappa
print '========'
print 'Mean P: ', meanP
print 'Error P:', errorP
print '========'
print 'Mean S: ', meanS/8.0
print 'Error S:', errorS/8.0
print '========'
print 'Mean U: ', (meanU-meanP)/8.0
print 'Error U:', (errorU+errorP)/8.0
print 'Mean U: ', (meanU)/8.0
print 'Error U:', (errorU)/8.0

## quit
## Emean = sum(E)/len(E)
## print Emean


            
                


