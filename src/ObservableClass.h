#ifndef OBSERVABLE_CLASS_H
#define OBSERVABLE_CLASS_H

#include "Common.h"
#include "PathDataClass.h"

class OutputFileClass
{
public:
  int dummy;
  /// Current has nothing in it.
};


/// This is the parent class for all observables.  It contains
/// a pointer to PathData.
class ObservableClass 
{
public:
  /// A reference to the PathData I'm observing
  PathDataClass &PathData;

  /// Observe the state of the present path and add it to the
  /// running sum for averages.
  virtual void Accumulate() = 0;
  /// Initialize my data.
  virtual void Initialize() = 0;
  /// Print running average to screen for debugging purposes
  virtual void Print() = 0;
  virtual void Write(OutputFileClass &outputFile)=0;
  /// The constructor.  Sets PathData references and calls initialize.
  ObservableClass(PathDataClass &myPathData) : PathData(myPathData)
  {   }
};



/// This template class will be used to construct distributed versions
/// of many different types of observables:  scalar observables, dVec
/// observables, array observables, array of dVec observables, etc.
/// We will write one class functions which correctly manages
/// collecting observables from all processors with MPI.
template<class T> 
class DistributedObservableClass : public ObservableClass
{
  int dummy;

  DistributedObservableClass(PathDataClass &myPathData) : ObservableClass (myPathData)
  { /* Do nothing for now. */ }
};

class ScalarObservableClass : public ObservableClass
{
protected:
  virtual void MyBlockAverage(double &mean, double &error) = 0;
public:
  /// This routine will collect averages of
  void Write (OutputFileClass &outputFile);

  ScalarObservableClass (PathDataClass &myPathData) : ObservableClass(myPathData)
  { /* Do nothing for now. */ }
};


/// A pair correlation function observable.
class PairCorrelation : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<int,1> Histogram;
  /// Stores the total number of counts
  int TotalCounts;
public:
  /// The species between which I am calculating the pair correlation
  /// function.
  int Species1, Species2;
  /// This grid defines the bins.  Bin 0 is bounded by 0 on the 
  /// bottom and grid(0) on the top.
  LinearGrid grid;
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  PairCorrelation(PathDataClass &myPathData) : ObservableClass(myPathData)
  { /* Do nothing for now. */ }
  void Write (OutputFileClass &outputFile)
    {/*Currently doesn't do anything*/};
};



#endif
