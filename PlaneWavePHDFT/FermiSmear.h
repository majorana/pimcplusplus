#ifndef FERMI_SMEAR_H
#define FERMI_SMEAR_H

#include "../Blitz.h"

class FermiSmearClass
{
public:
  virtual double D(double E, double mu) = 0;
  virtual double S(double E, double mu) = 0;
};

class MethfesselPaxton : public FermiSmearClass
{
protected:
  int Order;
  double Width;
  Array<double,1> A;
  Array<double,1> H;
  double S0(double x);
public:
  void   SetOrder (int order);
  inline int    GetOrder ()             { return Order; }
  inline void   SetWidth (double width) { Width = width; }
  inline double GetWidth ()             { return Width; }
  double D(double E, double mu);
  double S(double E, double mu);
};

#endif
  
