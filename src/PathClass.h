#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include "Common.h"
#include "MirroredArrayClass.h"

/*! A primitive class holding the postions of a group of identicle
  particles.  In addition, it holds the time in which each particle
  and timeslice was moved. */



class PathClass
{
 private:
  
  MirroredArrayClass<dVec> Positions;
  MirroredArrayClass<int> TimeStamp;
  
 public:
  /// Holds the MC step number each particle/timeslice was last updated.
  void shiftData(int slicesToShift,CommunicatorClass &Communicator);
  inline void resize(int numPtcles,int numTimeSlices)
    {
      Positions.resize(numPtcles,numTimeSlices);
      TimeStamp.resize(numPtcles,numTimeSlices);
    }
  
  /// Operator to access by value.
  inline dVec operator()(int Ptcl, int TimeSlice) const
    {
      return Positions(Ptcl, TimeSlice);
    }
  
  /// Write access function.
  inline void SetPos (int Ptcl, int TimeSlice, const dVec &NewPos)
    {
      
      Positions.Set(Ptcl,TimeSlice,NewPos);
      TimeStamp.Set(Ptcl,TimeSlice,GetCurrentTimeStamp());
    }
  
  inline int GetTimeStamp (int Ptcl, int TimeSlice)
    {
      return TimeStamp(Ptcl, TimeSlice);
    }
  

  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  void AcceptCopy (int Ptcl, int StartSlice, int EndSlice)
    {
      Positions.AcceptCopy(Ptcl,StartSlice,EndSlice);
      TimeStamp.AcceptCopy(Ptcl,StartSlice,EndSlice);
    }
  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  void RejectCopy (int Ptcl, int StartSlice, int EndSlice)
    {
      Positions.RejectCopy(Ptcl,StartSlice,EndSlice);
      TimeStamp.RejectCopy(Ptcl,StartSlice,EndSlice);
      
    }
};


#endif
