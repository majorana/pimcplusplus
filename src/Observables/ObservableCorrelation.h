#ifndef OBSERVABLE_CORRELATION_H
#define OBSERVABLE_CORRELATION_H

#include "ObservableBase.h"

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
  void WriteInfo();
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



///Creates a histogram of the end to end distance between the head and
///tail of the open loop.
class nofrClass : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<double,1> Histogram;
  //  //  Array<double,3> Histogram3d;
  //  //  ObservableVecDouble3 HistSum3d;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
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
  ///Writes Info about the observable when writing for the first time
  void WriteInfo();
  void Read(IOSectionClass& IO);
  nofrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection)
    ////    HistSum3d("3dSum",  IOSection,myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
  nofrClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection)
    ///    HistSum3d("3dSum",  IOSection,myPathData.Path.Communicator)
  { Initialize(); }

};






#endif
