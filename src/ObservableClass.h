#ifndef OBSERVABLE_CLASS_H
#define OBSERVABLE_CLASS_H

#include "Common.h"
#include "PathDataClass.h"

class ObservableClass
{
public:
  PathDataClass *PathData;

  virtual void Accumulate() = 0;
  void Initialize();
  void Print();
};


class PairCorrelation : public ObservableClass
{
  int Species1, Species2;
  Array<int,1> Histogram;
  int TotalCounts;
public:
  LinearGrid grid;
  void Accumulate();
  void Initialize();
  void Print();
};



#endif
