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

  IOSectionClass &IOSection;  
  virtual void WriteBlock()=0;
  virtual void ShiftData(int numTimeSlices)=0;
  /// The constructor.  Sets PathData references and calls initialize.
  ObservableClass(PathDataClass &myPathData,IOSectionClass &ioSection) : PathData(myPathData), IOSection(ioSection)
  {   }
};


/* class PrintConfigClass : public ObservableClass */
/* { */
/*  public: */
/*   void Accumulate(){;} */
/*   void Initialize(){;} */
/*   PrintConfigClass(PathDataClass &myPathData) : ObservableClass(myPathData){;} */
/*   void Write(OutputFileClass &outputFile){;} */
/*   void Print(){ */
/*     cout<<"Here are the current configurations"; */
/*     for (int counter=0;counter<PathData.NumSpecies();counter++){ */
/*       for (int counter2=0;counter2<PathData(counter).NumParticles();counter2++){ */
/* 	for (int counter3=0;counter3<PathData.NumTimeSlices();counter3++){ */
/* 	  cout<<PathData(counter3,counter,counter2)(0)<<" "; */
/* 	  cout<<PathData(counter3,counter,counter2)(1)<<" "; */
/* 	  cout<<PathData(counter3,counter,counter2)(2)<<" "; */
/* 	  cout<<endl; */
/* 	} */
/* 	cout<<endl; */
/*       } */
/*       cout<<endl; */
/*       cout<<endl; */
/*     } */
/*   } */
/* }; */
     

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

class LinkObservableClass : public ObservableClass
{
protected:
  //  virtual void MyBlockAverage(double &mean, double &error) = 0;
public:
  Array<double,1> LinkArray;
  double Total;
  int NumSamples;
  void Accumulate()
  {
    for (int counter=0;counter<LinkArray.size();counter++){
      Total+=LinkArray(counter);
    }
    NumSamples++;
  } 
  virtual void Update(int startTimeSlice,int endTimeSlice)=0;
  virtual void UpdateAll()=0;
  /// Initialize my data.
  virtual void Initialize() = 0;
  IOSectionClass &IOSection;  
  virtual void WriteBlock()=0;
  void ShiftData(int numTimeSlices);
  virtual void WriteBlock()=0;
  ScalarObservableClass (PathDataClass &myPathData, IOSectionClass &IOSection) : ObservableClass(myPathData), ObservableClass(IOSection)
  { 
    LinkArray.resize(PathData.NumTimeSlices());
    NumSamples=0;
  }
  
};

class EnergyClass : public LinkObservableClass
{
  

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
