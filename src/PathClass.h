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
  
  MirroredArrayClass<dVec> Positions; ///< The particle positions
  /// Holds the MC step number each particle/timeslice was last updated.
  MirroredArrayClass<int> TimeStamp;  
  
 public:
  /// Uses a Communicator to shift the date between processors so that
  /// the timeslices that each processor holds shifts by slicesToShift.
  void ShiftData(int slicesToShift,CommunicatorClass &communicator);
  /// Changes number of particles and timeslices;
  inline void Resize(int numPtcles,int numTimeSlices)
    {
      Positions.Resize(numPtcles,numTimeSlices);
      TimeStamp.Resize(numPtcles,numTimeSlices);
    }
  /// Returns the number of particles stored
  inline int NumParticles()
  {
    return Positions.NumParticles();
  }
  /// Returns the number of time slices
  inline int NumTimeSlices()
  {
    return Positions.NumTimeSlices();
  }
  /// Operator to access by value.
  inline dVec operator()(int ptcl, int timeSlice) const
    {
      return Positions(ptcl, timeSlice);
    }
  inline void MoveJoin(MirroredArrayClass1D<int> &PermMatrix,int oldJoin,int newJoin)
    {
      Positions.MoveJoin(PermMatrix,oldJoin,newJoin);
    }
  /// Write access function.
  inline void SetPos (int ptcl, int timeSlice, const dVec &newPos)
    {
      
      Positions.Set(ptcl,timeSlice,newPos);
      TimeStamp.Set(ptcl,timeSlice,GetCurrentTimeStamp());
    }
  /// Returns the MC stepnum in which ptcl/timeslice was last
  /// updated.
  inline int GetTimeStamp (int ptcl, int timeSlice)
    {
      return TimeStamp(ptcl, timeSlice);
    }
  

  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  void AcceptCopy (int ptcl, int startSlice, int endSlice)
    {
      Positions.AcceptCopy(ptcl,startSlice,endSlice);
      TimeStamp.AcceptCopy(ptcl,startSlice,endSlice);
    }
  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  void RejectCopy (int ptcl, int startSlice, int endSlice)
    {
      Positions.RejectCopy(ptcl,startSlice,endSlice);
      TimeStamp.RejectCopy(ptcl,startSlice,endSlice);
      
    }
};


#endif
