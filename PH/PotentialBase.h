#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#include "../IO/InputOutput.h"

class Potential
{
public:
  // Optional member functions -- if you're not a pseudoHamiltonian,
  // you do not need to define these
  virtual bool IsPH() { return false; }
  virtual double CoreRadius() { return 0.0; }
  virtual double A      (double r) { return 1.0; }
  virtual double B      (double r) { return 1.0; }
  virtual double dAdr   (double r) { return 0.0; }
  virtual double d2Adr2 (double r) { return 0.0; }

  // Required member functions
  virtual double V(double r) = 0;
  virtual double dVdr(double r) = 0;
  virtual double d2Vdr2(double r) = 0;
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
};


#endif
