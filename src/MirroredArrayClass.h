#ifndef MIRRORED_ARRAY_CLASS_H
#define MIRRORED_ARRAY_CLASS_H

#include "CommunicatorClass.h"

  
/// This is our magic backup array class.  It stores two copies of a
///one-dimensional array in a two-dimensional internal array, an
///active copy and a backed-up copy.  The first copy is indicated by a
///first index of 0 and the second by 1.  It has two modes, MOVE and
///OBSERVABLE.  In MOVE mode, it writes only to the first copy, while
///in observable mode, it writes to both copies.  The AcceptCopy and
///RejectCopy member functions are used to accept or reject moves.
///Currently, I have pulled out the shiftData, because it wasn't exactly clear 
///to me what exactly it meant to shift the data in this case

template<class T>
class MirroredArrayClass1D
{
private:

  /// Array holds the A and B copies of a two dimensional array 
  Array<T,2> AB; /// (0=A 1=B, particles)



public:
  /// Resizes the two dimensional array.
  inline void Resize(int numPtcls)
  {  AB.resize(2,numPtcls); }

  inline int NumParticles()
  {
    return AB.extent(1);
  }

  /// Constructor that creates the 2d array of the correct size
  inline MirroredArrayClass1D(int particleNum)
  {
    AB.resize(2,particleNum);
  }
  /// Constructor that does nothing.
  MirroredArrayClass1D(){};
  /// Debug Printing of some sort
  void  Print();
  /// Returns the active value. 
  inline T operator()(int x) const
  {
    return AB(Write1,x);
  }


  
  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int x, const T &newVal)
  {
    AB(Write1,x)=newVal;
    AB(Write2,x)=newVal;
  }


  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  inline void AcceptCopy (int particle)
  {
      AB(1,particle) = AB(0,particle);
  }

  inline void AcceptCopy (const Array<int,1> &activePtcls)
  {
    for (int i=0; i<activePtcls.size(); i++)
      AcceptCopy(activePtcls(i));
  }

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  inline void RejectCopy (int particle)
  {
      AB(0,particle) = AB(1,particle);
  }
  inline void RejectCopy (const Array<int,1> &activePtcls)
  {
    for (int i=0; i<activePtcls.size(); i++)
      RejectCopy(activePtcls(i));
  }


};

 
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
  Array<T,3> AB; /// (0=A 1=B, timeslice,particles)

public:
  /// Resizes the two dimensional array.
  inline void Resize(int numTimeSlices,int numPtcles)
  { AB.resize(2,numTimeSlices,numPtcles); } 
  void MoveJoin(MirroredArrayClass1D<int> &PermMatrix,
		int oldJoin, int newJoin);
  inline int NumParticles()
  {
    return AB.extent(2);
  }
  inline int NumTimeSlices()
  {
    return AB.extent(1);
  }
  /// Constructor that creates the 2d array of the correct size
  inline MirroredArrayClass(int timeSliceNum, int particleNum)
  {
    AB.resize(2,timeSliceNum,particleNum);
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
  void ShiftData(int slicesToShift,PIMCCommunicatorClass &communicator); 

  
  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int x, int y,const T &newVal)
  {
    AB(Write1,x,y)=newVal;
    AB(Write2,x,y)=newVal;
  }

  inline void AcceptCopy (int slice, int ptcl){
    AB(1,slice,ptcl)=AB(0,slice,ptcl);
  }
  inline void RejectCopy(int slice,int ptcl){
    AB(0,slice,ptcl)=AB(1,slice,ptcl);
  }
    


  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  inline void AcceptCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++){
      for (int i=0; i<activeParticles.size(); i++){
	int ptcl=activeParticles(i);
	AB(1,slice,ptcl) = AB(0,slice,ptcl);	
      }
    }
  }

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  inline void RejectCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++){
      for (int i=0; i<activeParticles.size(); i++){
	int ptcl=activeParticles(i);
	AB(0,slice,ptcl) = AB(1,slice,ptcl);	
      }
    }
  }

 
};


/// This class holds a mirrored tensor of rank (TimeSlice x Ptcl x
/// Ptcl).  The (Ptcl x Ptcl) matrix is symmetric, so that we store
/// only the lower diagonal part of the matrix (including the diagonal
/// terms.)  It can be accessed with functions that take 
/// (timeSlice, ptcl1, ptcl2) or (timeSlice, PairIndex).

template<class T>
class MirroredSymmetricMatrixClass
{
private:
  /// Array holds the A and B copies of a two dimensional array 
  Array<T,3> AB; ///< (0=A 1=B, timeslice, pairIndex)
  /// This is the index into the array.  We map (ptcl #)x(ptcl #)
  /// into a single integer
  int CurrentPairNum;
  int NumPtcls;
  inline void Order (int &i1, int &i2) const
  {
    if (i2 > i1)
      {
	int temp = i2;
	i2 = i1;
	i1 = temp;
      }
  }

  inline int PairIndex(int ptcl1, int ptcl2) const
    { 
      Order(ptcl1, ptcl2);
      return (((ptcl1*(ptcl1+1))>>1)+ptcl2); 
    }
public:
  /// Resizes the two dimensional array.
  inline void Resize(int numTimeSlices,int numPtcls)
  {
    NumPtcls = numPtcls;
    int NumPairs = ((numPtcls*(numPtcls+1))>>1);
    AB.resize(2, numTimeSlices, NumPairs);
  }

  inline int NumParticles()
  {
    return NumPtcls;
  }
  inline int NumTimeSlices()
  {
    return AB.extent(1);
  }
  /// Debug Printing of some sort
  void  Print();
  /// Returns the active value. 
  inline T operator()(int timeSlice, int ptcl1, int ptcl2) const
  {
    Order(ptcl1, ptcl2);
    int index = PairIndex(ptcl1, ptcl2);
    return AB(Write1,timeSlice,index);
  }

  inline T operator()(int timeSlice, int pairIndex) const
  {
    return AB(Write1, timeSlice, pairIndex);
  }

  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int timeSlice, int pairIndex ,const T &newVal)
  {
    AB(Write1,timeSlice,pairIndex)=newVal;
    AB(Write2,timeSlice,pairIndex)=newVal;
  }

  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int timeSlice, int ptcl1, int ptcl2 ,const T &newVal)
  {
    Order(ptcl1, ptcl2);
    int index = PairIndex(ptcl1,ptcl2);
    AB(Write1,timeSlice,index)=newVal;
    AB(Write2,timeSlice,index)=newVal;
  }

  ///This shifts slicesToShift time slices to the next (if positive)
  ///or previous (if negative) processor 
  void ShiftData(int slicesToShift,PIMCCommunicatorClass &communicator); 

  void MoveJoin(MirroredArrayClass1D<int> &PermMatrix,
		int oldJoin, int newJoin);
      
  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  inline void AcceptCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++) {
      for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
	int ptcl1=activeParticles(ptcl1Index);
	for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++) {
	  int index = PairIndex(ptcl1,ptcl2);
	  AB(1,slice,index) = AB(0,slice,index);	
	}
      }
    }
  }    

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  inline void RejectCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++) {
      for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
	int ptcl1=activeParticles(ptcl1Index);
	for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++) {
	  int index = PairIndex(ptcl1,ptcl2);
	  AB(0,slice,index) = AB(1,slice,index);	
	}
      }
    }
  }    

 
  /// Constructor that creates the 2d array of the correct size
  inline MirroredSymmetricMatrixClass(int numTimeSlices, int numPtcls)
  {
    Resize(numTimeSlices, numPtcls);
  }
  /// Constructor that does nothing.
  MirroredSymmetricMatrixClass(){};


};






/// This class holds a mirrored tensor of rank (TimeSlice x Ptcl x
/// Ptcl).  The (Ptcl x Ptcl) matrix is antisymmetric, so that we store
/// only the lower diagonal part of the matrix (including the diagonal
/// terms.)  It can be accessed with functions that take 
/// (timeSlice, ptcl1, ptcl2) or (timeSlice, PairIndex).
/// Automatically reverses the sign if the ptcl1 < ptcl2.

template<class T>
class MirroredAntiSymmetricMatrixClass
{
private:
  /// Array holds the A and B copies of a two dimensional array 
  Array<T,3> AB; ///< (0=A 1=B, timeslice, pairIndex)
  /// This is the index into the array.  We map (ptcl #)x(ptcl #)
  /// into a single integer
  int CurrentPairNum;
  int NumPtcls;
  
  inline void swap (int &ptcl1, int &ptcl2) const
  {
    int temp = ptcl1;
    ptcl1 = ptcl2;
    ptcl2 = temp;
  }

  inline int PairIndex(int ptcl1, int ptcl2) const
  { 
    if (ptcl2 > ptcl1)
      swap(ptcl1,ptcl2);
    return (((ptcl1*(ptcl1+1))>>1)+ptcl2); 
  }
  inline int PairIndex(int ptcl1, int ptcl2, bool swapSign) const
  { 
    if (ptcl2 > ptcl1) {
      swap(ptcl1,ptcl2);
      swapSign = true;
    }
    else
      swapSign=false;     
    return (((ptcl1*(ptcl1+1))>>1)+ptcl2); 
  }
public:
  /// Resizes the two dimensional array.
  inline void Resize(int numTimeSlices,int numPtcls)
  {
    NumPtcls = numPtcls;
    int NumPairs = ((numPtcls*(numPtcls+1))>>1);
    AB.resize(2, numTimeSlices, NumPairs);
  }

  inline int NumParticles()
  {
    return NumPtcls;
  }
  inline int NumTimeSlices()
  {
    return AB.extent(1);
  }
  /// Debug Printing of some sort
  void  Print();
  /// Returns the active value. 
  inline T operator()(int timeSlice, int ptcl1, int ptcl2) const
  {
    bool swapSign;
    int index = PairIndex(ptcl1, ptcl2, swapSign);
    if (swapSign)
      return (-AB(Write1, timeSlice, index));
    else 
      return (AB(Write1, timeSlice, index));
  }

  inline T operator()(int timeSlice, int pairIndex) const
  {
    return AB(Write1, timeSlice, pairIndex);
  }

  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int timeSlice, int pairIndex ,const T &newVal)
  {
    AB(Write1,timeSlice,pairIndex)=newVal;
    AB(Write2,timeSlice,pairIndex)=newVal;
  }

  /// Write to the mirrored array in the way specified by the present
  /// mode. 
  inline void Set(int timeSlice, int ptcl1, int ptcl2 ,const T &newVal)
  {
    if (ptcl2 > ptcl1)
      {
	swap (ptcl1, ptcl2);
	int index = PairIndex(ptcl1, ptcl2);
	AB(Write1, timeSlice, index) = -newVal;
	AB(Write2, timeSlice, index) = -newVal;
      }
    else
      {
	int index = PairIndex(ptcl1, ptcl2);
	AB(Write1,timeSlice,index) = newVal;
	AB(Write2,timeSlice,index) = newVal;
      }
  }

  ///This shifts slicesToShift time slices to the next (if positive)
  ///or previous (if negative) processor 
  void ShiftData(int slicesToShift,PIMCCommunicatorClass &communicator); 

  void MoveJoin(MirroredArrayClass1D<int> &PermMatrix,
		int oldJoin, int newJoin);
      
  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from A to B.
  inline void AcceptCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++) {
      for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
	int ptcl1=activeParticles(ptcl1Index);
	for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++) {
	  int index = PairIndex(ptcl1,ptcl2);
	  AB(1,slice,index) = AB(0,slice,index);	
	}
      }
    }
  }    

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  inline void RejectCopy (int startSlice, int endSlice, 
			  const Array<int,1> &activeParticles)
  {
    for (int slice=startSlice; slice<=endSlice; slice++) {
      for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
	int ptcl1=activeParticles(ptcl1Index);
	for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++) {
	  int index = PairIndex(ptcl1,ptcl2);
	  AB(0,slice,index) = AB(1,slice,index);	
	}
      }
    }
  }    

 
  /// Constructor that creates the 2d array of the correct size
  inline MirroredAntiSymmetricMatrixClass(int numTimeSlices, int numPtcls)
  {
    Resize(numTimeSlices, numPtcls);
  }
  /// Constructor that does nothing.
  MirroredAntiSymmetricMatrixClass(){};


};



#endif
