#ifndef BACKUP_ARRAY_CLASS_H
#define BACKUP_ARRAY_CLASS_H

/*! This is our magic backup array class.  It stores two copies of a
two-dimensional array in a three-dimensional internal array, an active
copy and a backed-up copy.  The first copy is indicated by a first
index of 0 and the second by 1.  It has two modes, MOVE and
OBSERVABLE.  In MOVE mode, it writes only to the first copy, while in
observable mode, it writes to both copies.  The AcceptCopy and
RejectCopy member functions are used to accept or reject moves. */
class MirroredArrayClass<T>
{
private:
  Array<T,3> AB;
  int Write1, Write2;

  /// Returns the active value.
  inline T operator()(int x,int y) const
  {
    return AB(0,x,y);
  }

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

  /// In MOVE mode, write only to the first copy, while in OBSERVABLE
  /// mode, write to both copies.
  inline SetMode(ModeType Mode)
  {
    /// Write only to the first copy
    if (ModeType == MOVEMODE) 
      {
	Write1 = 0;
	Write2 = 0;
      }
    else if (ModeType == OBSERVABLEMODE)
      {
	Write1 = 0;
	Write2 = 1;
      }
    else
      {
	abort ("Undefined mode type!");
      }
  }

  /// In case of acceptance, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from 0 to 1.
  void AcceptCopy (const Array<int,1> &Ptcls, int StartSlice, int EndSlice)
  {
    for (int i=0; i<Ptcls.size(); i++)
      for (int j=StartSlice; j<=EndSlice; j++)
	AB(1,i,j) = AB(0,i,j);
  }

  /// In case of rejection, this is called to copy the new path over
  /// the backup copy.  StartSlice and EndSlice are inclusive.  This
  /// copies from 1 to 0.
  void RejectCopy (const Array<int,1> &Ptcls, int StartSlice, int EndSlice)
  {
    for (int i=0; i<Ptcls.size(); i++)
      for (int j=StartSlice; j<=EndSlice; j++)
	  AB(0,i,j) = AB(1,i,j);
  }


};
