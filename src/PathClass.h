#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include "Blitz.h"
/*! A primitive class holding the postions of a group of identicle
  particles.  In addition, it holds the time in which each particle
  and timeslice was moved. */
class PathClass
{
private:
  Array<dVec, 2> Positions;

public:
  /// Holds the MC step number each particle/timeslice was last updated.
  Array<int, 2> TimeStamp;

  /// Operator to access by value.
  inline dVec operator()(int Ptcl, int TimeSlice) const
  {
    return Positions(Ptcl, TimeSlice);
  }
  /// Operator to access by reference.
  inline dVec& operator()(int Ptcl, int TimeSlice)
  {
    return Positions(Ptcl, TimeSlice);
  }
}


#endif
