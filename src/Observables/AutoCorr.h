#ifndef AUTO_CORRELATION_H
#define AUTO_CORRELATION_H

#include "ObservableBase.h"

/// An autocorrelation observable.
class AutoCorrClass : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<double,1> Histogram;
  ObservableVecDouble1 HistVar;
  /// Stores dipole moments
  Array<dVec,2> DipoleBin;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  int now;
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
  void WriteInfo();
  void Read(IOSectionClass& IO);
  int MapIndex(int slice, int molecule);
  void Advance(int& index,int limit);
  dVec MeasureDipole(int slice,int molecule);
  double CalcAutoCorr(int index, int t, int limit);
  int Locate(int i, int t, int limit);
  double CalcDotProd(dVec v1, dVec v2);
  dVec Rotate(dVec coord, double theta);
  AutoCorrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    HistVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
    Initialize();
  }

/*  AutoCorrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData, ioSection) 
  { Initialize(); }
*/

};
#endif
