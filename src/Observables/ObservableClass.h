#ifndef OBSERVABLE_CLASS_H
#define OBSERVABLE_CLASS_H

#include "ObservableEnergy.h"
#include "ObservableCorrelation.h"
#include "PathDump.h"
#include "WindingNumber.h"
#include "StructureFactor.h"     
#include "Weight.h"
#include "ObservableModifiedEnergy.h"
#include "DistanceToHead.h"
#include "Time.h"
#include "Coupling.h"
#include "VacancyLocation.h"
#include "AutoCorr.h"
#include "Angular.h"
#include "AngularMomCor.h"
#include "SuperfluiDrop.h"
#include "VacancyLocation.h"
#include "Forces.h"
#include "Phik.h"
#include "JosephsonPathDump.h"
#include "PermutationCount.h"
/// This template class will be used to construct distributed versions
/// of many different types of observables:  scalar observables, dVec
/// observables, array observables, array of dVec observables, etc.
/// We will write one class functions which correctly manages
/// collecting observables from all processors with MPI.
// template<class T> 
// class DistributedObservableClass : public ObservableClass
// {
//   int dummy;

//   DistributedObservableClass(PathDataClass &myPathData) : ObservableClass (myPathData)
//   { /* Do nothing for now. */ }
// };








#endif
