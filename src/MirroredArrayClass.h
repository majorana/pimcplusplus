#ifndef MIRRORED_ARRAY_CLASS_H
#define MIRRORED_ARRAY_CLASS_H

#include "CommunicatorClass.h"

/*! This is our magic backup array class.  It stores two copies of a
two-dimensional array in a three-dimensional internal array, an active
copy and a backed-up copy.  The first copy is indicated by a first
index of 0 and the second by 1.  It has two modes, MOVE and
OBSERVABLE.  In MOVE mode, it writes only to the first copy, while in
observable mode, it writes to both copies.  The AcceptCopy and
RejectCopy member functions are used to accept or reject moves. */
template<class T>
class MirroredArrayClass
{
private:
  Array<T,3> AB; /// (0=A 1=B, particles, timeslice)



public:
  void resize(int numPtcles,int numTimeSlices);
  inline MirroredArrayClass(int particleNum, int timeSliceNum)
  {
    AB.resize(2,particleNum,timeSliceNum);
  }
  MirroredArrayClass(){};
  void  Print();
  /// Returns the active value.
  inline T operator()(int x,int y) const
  {
    return AB(Write1,x,y);
  }
  ///This shifts slicesToShift time slices to the next (if positive) or previous (if negative) processor
  void shiftData(int slicesToShift,CommunicatorClass &Communicator); 
  /// Returns the backup value.
  inline T Backup (int x, int y) const
  {
    return (1,x,y);
  }
  
  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int x, int y,const T &NewVal)
  {
    AB(Write1,x,y)=NewVal;
    AB(Write2,x,y)=NewVal;
  }


  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from 0 to 1.
  inline void AcceptCopy (int Particle, int StartSlice, int EndSlice)
  {
    for (int Slice=StartSlice; Slice<=EndSlice; Slice++)
      AB(1,Particle,Slice) = AB(0,Particle,Slice);
  }

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from 1 to 0.
  inline void RejectCopy (int Particle, int StartSlice, int EndSlice)
  {
    for (int Slice=StartSlice; Slice<=EndSlice; Slice++)
      AB(0,Particle,Slice) = AB(1,Particle,Slice);
  }
};




#endif
