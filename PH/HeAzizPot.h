#ifndef HEAZIZ_POT_H
#define HEAZIZ_POT_H

#include "PotentialBase.h"

class HeAzizPot : public Potential
{
private:
  static const double bb[7];
public:
  
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);
};

#endif
