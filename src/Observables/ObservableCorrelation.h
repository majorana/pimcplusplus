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



/// A pair correlation function observable.
class nofrClass : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<double,1> Histogram;
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
  nofrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection){
    TimesCalled=0;
  }
  nofrClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection), 
    Species1(species1), Species2(species2) 
  { Initialize(); }

};

///Need to be moved out in a minute or two
class WindingNumberClass : public ObservableClass
{
 private:
  int Freq;
  int DumpFreq;
  dVec TotalW2;
  Array<dVec,1> TotalDisp;
  Array<dVec,1> TempDisp;
 public:
  int TimesCalled;
  int NumSamples;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  WindingNumberClass(PathDataClass &myPathData,IOSectionClass &ioSection):
    ObservableClass(myPathData,ioSection) {
    TimesCalled=0;
    NumSamples=0;
    TotalDisp.resize(PathData.Path.NumParticles());
    TempDisp.resize(PathData.Path.NumParticles());
    TotalW2=0.0;
  }




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
