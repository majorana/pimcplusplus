#ifndef OBSERVABLE_CLASS_H
#define OBSERVABLE_CLASS_H

#include "Common.h"
#include "PathDataClass.h"
#include "Common/IO/InputOutput.h"


/// This is the parent class for all observables.  It contains
/// a pointer to PathData.
class ObservableClass 
{
protected:
  bool FirstTime;
  VarClass *IOVar;
public:
  /// A reference to the PathData I'm observing
  PathDataClass &PathData;
  /// Note: This is not a reference.  If it were, it could change
  /// behind our backs
  IOSectionClass IOSection;  
  string Name;
  /// Observe the state of the present path and add it to the
  /// running sum for averages.
  virtual void Accumulate() = 0;
  virtual void WriteBlock()=0;
  virtual void Read(IOSectionClass& IO)=0;
  virtual void ShiftData(int numTimeSlices) {;}
  /// The constructor.  Sets PathData references and calls initialize.
  ObservableClass(PathDataClass &myPathData,IOSectionClass ioSection) 
    : PathData(myPathData), IOSection(ioSection)
  {
    FirstTime = true;
    Name="";
  }
};


     

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




/// A pair correlation function observable.
class PairCorrelationClass : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<int,1> Histogram;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
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
  void WriteBlock();
  void Read(IOSectionClass& IO);
  PairCorrelationClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection){
    TimesCalled=0;
  }
  PairCorrelationClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection), 
    Species1(species1), Species2(species2) 
  { Initialize(); }


};


class PathDumpClass : public ObservableClass
{
private:
public:
  int TimesCalled;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  PathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection)  { 
    Name="PathDump";
    TimesCalled=0;
  }
};




#endif
