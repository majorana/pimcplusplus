#ifndef BACKUP_ARRAY_CLASS_H
#define BACKUP_ARRAY_CLASS_H

/*! This is our magic backup array class */

class MirroredArrayClass<T>
{
private:
  Array<T,2> A;
  Array<T,2> B;
  Array<T,2> *Write1, *Write2;

  inline T operator()(int x,int y) const
  {
    return A(x,y);
  }

  inline void Set(int x, int y,const T &NewVal)
  {
    (*Write1)(x,y)=NewVal;
    (*Write2)(x,y)=NewVal;
  }

  

  inline SetMode(ModeType Mode)
  {
    if (ModeType == MOVEMODE) 
      {
	Write1 = &A;
	Write2 = &A;
      }
    else if (ModeType == OBSERVABLEMODE)
      {
	Write1 = &A;
	Write2 = &B;
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
	B(i,j) = A(i,j);

       
  }
  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from B to A.
  void RejectCopy (const Array<int,1> &Ptcls, int StartSlice, int EndSlice)
  {
    for (int i=0; i<Ptcls.size(); i++)
      for (int j=StartSlice; j<=EndSlice; j++)
	  A(i,j) = B(i,j);

      
  }


};
