#ifndef SUPERFLUID_FRACTION__H
#define SUPERFLUID_FRACTION__H

#include "ObservableBase.h"
#include "WindingNumber.h"



class SuperfluidFractionClass : public WindingNumberClass
{
 protected:
  ObservableVecDouble1 SFVar;
 public:
  void Read(IOSectionClass& IO);
  void WriteBlock();
  SuperfluidFractionClass(PathDataClass &myPathData,IOSectionClass &ioSection) :
    WindingNumberClass(myPathData,ioSection),
    SFVar("SuperfluidFraction", IOSection, myPathData.Path.Communicator)
  {
    W2Sum = 0.0;
    SamplesInBlock = 0;
  }

};



#endif
