#!/usr/bin/python
import tables
import sys
import stats
#import Numeric
h5file = tables.openFile(sys.argv[1])
#data = h5file.root.Paths_2.Path.read()
#names=h5file.root.Paths_2.SpeciesNames.read()

E = h5file.root.Energies_0.TotalEnergy.read()
(mean,var,error,kappa)= stats.Stats(E)
print 'Mean:  ', mean
print 'Var:   ', var
print 'Error: ', error
print 'kappa  ', kappa
## quit
## Emean = sum(E)/len(E)
## print Emean


            
                


