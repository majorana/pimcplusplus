#ifndef MIRRORED_ARRAY_CLASS_H
#define MIRRORED_ARRAY_CLASS_H

#include "CommunicatorClass.h"

/// This is our magic backup array class.  It stores two copies of a
///two-dimensional array in a three-dimensional internal array, an
///active copy and a backed-up copy.  The first copy is indicated by a
///first index of 0 and the second by 1.  It has two modes, MOVE and
///OBSERVABLE.  In MOVE mode, it writes only to the first copy, while
///in observable mode, it writes to both copies.  The AcceptCopy and
///RejectCopy member functions are used to accept or reject moves.

template<class T>
class MirroredArrayClass
{
private:

  /// Array holds the A and B copies of a two dimensional array 
  Array<T,3> AB; /// (0=A 1=B, particles, timeslice)



public:
  /// Resizes the two dimensional array.
  void Resize(int numPtcles,int numTimeSlices);

  inline int NumParticles()
  {
    return AB.extent(1);
  }
  inline int NumTimeSlices()
  {
    return AB.extent(2);
  }
  /// Constructor that creates the 2d array of the correct size
  inline MirroredArrayClass(int particleNum, int timeSliceNum)
  {
    AB.resize(2,particleNum,timeSliceNum);
  }
  /// Constructor that does nothing.
  MirroredArrayClass(){};
  /// Debug Printing of some sort
  void  Print();
  /// Returns the active value. 
  inline T operator()(int x,int y) const
  {
    return AB(Write1,x,y);
  }

  ///This shifts slicesToShift time slices to the next (if positive)
  ///or previous (if negative) processor 
  void ShiftData(int slicesToShift,CommunicatorClass &communicator); 

  
  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int x, int y,const T &newVal)
  {
    AB(Write1,x,y)=newVal;
    AB(Write2,x,y)=newVal;
  }


  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  inline void AcceptCopy (int particle, int startSlice, int endSlice)
  {
    for (int slice=startSlice; slice<=endSlice; slice++)
      AB(1,particle,slice) = AB(0,particle,slice);
  }

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  inline void RejectCopy (int particle, int startSlice, int endSlice)
  {
    for (int slice=startSlice; slice<=endSlice; slice++)
      AB(0,particle,slice) = AB(1,particle,slice);
  }
};




#endif
