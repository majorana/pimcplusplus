#ifndef BCCLASS_H
#define BCCLASS_H
#include "Common.h"

class BCClass 
{
private:
  PathClass &Path;
  dVec IsPeriodic;
  dVec Box, BoxInv;
public:
  inline void SetBox (dVec box);
  inline dVec GetBox ();
  
  inline void SetPeriodic(TinyVector<bool,NDIM> period);

  inline void DistDisp (int slice, int ptcl1, int ptcl2,
			double &dist, dVec &disp);

  inline void DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
			double &distA, double &distB,dVec &dispA, dVec &dispB);

  inline double Distance (int slice, int ptcl1, int ptcl2);

  inline dVec Velocity (int sliceA, int sliceB, int ptcl);

  BCClass (PathClass &path) : Path(path)
  { /* do nothing */ }
};

inline void BCClass::SetBox (dVec box)
{
  Box = box;
  for (int i=0; i<NDIM; i++)
    BoxInv(i) = 1.0/box(i);
}

inline dVec BCClass::GetBox ()
{
  return Box;
}


inline void BCClass::SetPeriodic(TinyVector<bool,NDIM> period)
{
  for (int i=0; i<NDIM; i++)
    IsPeriodic(i) = period(i) ? 1.0 : 0.0;
}

inline void BCClass::DistDisp (int slice, int ptcl1, int ptcl2,
			       double &dist, dVec &disp)
{
  disp = Path(slice, ptcl2) -Path(slice, ptcl1);
  
  for (int i=0; i<NDIM; i++) {
    double n = -trunc(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl2) -Path(slice, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i))
      DBdisp(i) += Box(i);
    assert (fabs(DBdisp(i)-disp(i)) < 1.0e-12);
  }
#endif
}


inline void BCClass::DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
			       double &distA, double &distB, 
			       dVec &dispA, dVec &dispB)
{  
  dispA = Path(sliceA, ptcl2) - Path(sliceA,ptcl1);
  dispB = Path(sliceB, ptcl2) - Path(sliceB,ptcl1);
  
  for (int i=0; i<NDIM; i++) {
    double n = -trunc(dispA(i)*BoxInv(i)+0.5);
    dispA(i) += n*IsPeriodic(i)*Box(i);
    dispB(i) += n*IsPeriodic(i)*Box(i);
  }
  distA = sqrt(dot(dispA,dispA));
  distB = sqrt(dot(dispB,dispB));

#ifdef DEBUG
  dVec DBdispA = Path(sliceA, ptcl2) -Path(sliceA, ptcl1);
  dVec DBdispB = Path(sliceB, ptcl2) -Path(sliceB, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdispA(i) > 0.5*Box(i)) {
      DBdispA(i) -= Box(i);
      DBdispB(i) -= Box(i);
    }
    while (DBdispA(i) < -0.5*Box(i)) {
      DBdispA(i) += Box(i);
      DBdispB(i) += Box(i);
    }
    assert (fabs(DBdispA(i)-dispA(i)) < 1.0e-12);
    assert (fabs(DBdispB(i)-dispB(i)) < 1.0e-12);
  }
#endif
}

inline dVec BCClass::Velocity (int sliceA, int sliceB, int ptcl)
{
  return Path(sliceB, ptcl) - Path(sliceA,ptcl);
}


#endif
