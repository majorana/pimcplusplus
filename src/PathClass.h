#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include "Common.h"

/*! A primitive class holding the postions of a group of identicle
  particles.  In addition, it holds the time in which each particle
  and timeslice was moved. */
class PathClass
{
private:
  Array<dVec, 2> PositionsA, PositionsB;
  Array<int, 2> TimeStampA, TimeStampB;
  Array<int, 2> *WriteTimeStamp1, *WriteTimeStamp2;
  Array<dVec, 2> *Write1, *Write2;
public:
  /// Holds the MC step number each particle/timeslice was last updated.


  /// Operator to access by value.
  inline dVec operator()(int Ptcl, int TimeSlice) const
  {
    return PositionsA(Ptcl, TimeSlice);
  }

  /// Write access function.
  inline void SetPos (int Ptcl, int TimeSlice, const dVec &NewPos)
  {
    (*Write1)(Ptcl, TimeSlice) = NewPos;
    (*Write2)(Ptcl, TimeSlice) = NewPos;
    (*WriteTimeStamp1)(Ptcl, TimeSlice) = GetCurrentTimeStamp();
    (*WriteTimeStamp2)(Ptcl, TimeSlice) = GetCurrentTimeStamp();
  }

  inline TimeStamp (int Ptcl, int TimeSlice)
  {
    return TimeStampA(Ptcl, TimeSlice);
  }

  inline SetMode(ModeType Mode)
  {
    if (ModeType == MOVEMODE) 
      {
	Write1 = &PositionsA;
	Write2 = &PositionsA;
      }
    else if (ModeType == OBSERVABLEMODE)
      {
	Write1 = &PositionsA;
	Write2 = &PositionsB;
      }
    else
      {
	abort ("Undefined mode type!");
      }
  }
  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  void AcceptCopy (const Array<int,1> &Ptcls, int StartSlice, int EndSlice)
  {
    for (int i=0; i<Ptcls.size(); i++)
      for (int j=StartSlice; j<=EndSlice; j++)
	{
	  PositionsB(i,j) = PositionsA(i,j);
	  TimeStampB(i,j) = TimeStampA(i,j);
	}
  }
  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  void RejectCopy (const Array<int,1> &Ptcls, int StartSlice, int EndSlice)
  {
    for (int i=0; i<Ptcls.size(); i++)
      for (int j=StartSlice; j<=EndSlice; j++)
	{
	  PositionsA(i,j) = PositionsB(i,j);
	  TimeStampA(i,j) = TimeStampB(i,j);
	}
  }
}


#endif
