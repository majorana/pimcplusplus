#ifndef QUINTIC_PH_H
#define QUINTIC_PH_H

#include "PotentialBase.h"
#include "../Splines/QuinticSpline.h"

class QuinticPH : public Potential
{
private:
  QuinticSpline pA, pB, Vcore;
  Potential *Vouter;
  double CoreRadius;
  LinearGrid Agrid, Bgrid, Vgrid;
  // The minimum value of A and B
  double ABmin;
public:

  // This is true if this represents a bare potential.  Otherwise, it
  // represents a total potential:  Total = Bare + XC
  bool IsBare;
  // True if we want to use Vcore spline
  bool UseVcore;

  bool IsPH() { return true; }

  inline double GetCoreRadius ();
  inline void   SetCoreRadius (double coreRadius);

  // These are used to get and set the minimum value of A and B,
  // adjusting PA and PB accordingly.
  inline double GetABmin();
  inline double SetABmin();

  // These are used to get and set the values of A, B and V
  inline int    NumAParams();
  inline int    NumBParams();
  inline int    NumVParams();
  inline double GetAval (int i);
  inline void   SetAval (int i, double Aval, bool isNegative);
  inline double GetBval (int i);
  inline void   SetBval (int i, double Bval, bool isNegative);
  inline double GetVval (int i);
  inline void   SetVval (int i, double val);

  double A      (double r);
  double dAdr   (double r);
  double d2Adr2 (double r);
  double B      (double r);
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  inline double  GetParam (int i) const;
  inline void SetParam (int i, double val);
  inline int NumParams();
  void Read (IOSectionClass &in);
  void Write (IOSectionClass &out); 
};

inline double QuinticPH::GetCoreRadius()
{
  return (CoreRadius);
}

inline void QuinticPH::SetCoreRadius(double coreRadius)
{
  CoreRadius = coreRadius;
}

inline double QuinticPH::GetParam(int i) const
{
  // A is constrained at core radius
  int numA = pA.NumPoints()-1;
  // B is constrained at origin at core radius
  int numB = pB.NumPoints()-2;
  // V is constrained at core radius
  int numV = Vcore.NumPoints()-1;
  if (i < numA)
    return (pA(i));
  else if (i < (numA+numB))
    return (pB(i-numA+1));
  else 
    return (Vcore(i-numA-numB));
}

inline void QuinticPH::SetParam(int i, double val)
{
  // A is constrained at core radius
  int numA = pA.NumPoints()-1;
  // B is constrained at origin at core radius
  int numB = pB.NumPoints()-2;
  // V is constrained at core radius
  int numV = Vcore.NumPoints()-1;

  // B(0) is locked to be A(0)
  if (i == 0) {
    pA(0) = val; 
    pB(0) = val;
  }
  else if (i < numA)
    pA(i) = val;
  else if (i < (numA+numB))
    pB(i-numA+1) = val;
  else 
    Vcore(i-numA-numB) = val;
}

    
inline double QuinticPH::GetAval (int i)
{
  return (pA(i)*pA(i) + ABmin);
}

inline void QuinticPH::SetAval (int i, double Aval, bool isNegative)
{
  double pAval = sqrt (Aval - ABmin);
  if (isNegative)
    pAval = -pAval;

  pA(i) = Aval;
}


inline double QuinticPH::GetBval (int i)
{
  return (pA(i)*pA(i) + ABmin);
}

inline void QuinticPH::SetBval (int i, double Bval, bool isNegative)
{
  double pBval = sqrt (Bval - ABmin);
  if (isNegative)
    pBval = -pBval;

  pB(i+1) = Bval;
}

inline double QuinticPH::GetVval (int i)
{
  return (Vcore(i));
}

inline void QuinticPH::SetVval (int i, double val)
{
  Vcore(i) = val;
}

inline int QuinticPH::NumAParams()
{
  return (pA.NumPoints()-1);
}

inline int QuinticPH::NumBParams()
{
  return (pB.NumPoints()-2);
}

inline int QuinticPH::NumVParams()
{
  return (Vcore.NumPoints()-1);
}


#endif
