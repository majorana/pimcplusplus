#!/usr/bin/python
from tables import *

h5file = openFile("../Observables2.h5");
data = h5file.root.Energies_0.TotalEnergy.read()

print data
