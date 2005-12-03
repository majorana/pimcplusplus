#ifndef MIRRORED_CLASS
#define MIRRORED_CLASS

#include <blitz/array.h>

typedef enum {OLDMODE, NEWMODE} ModeType;
extern ModeType ActiveCopy;

inline void SetMode (ModeType type)
{
  ActiveCopy=type;
  return;
  //  if (type == OLDMODE)
  //    ActiveCopy = 0;
  //  else 
  //    ActiveCopy = 1;
}

inline ModeType GetMode()
{
  return ActiveCopy;
}

using namespace blitz;

template <class T>
class MirroredClass
{
private:
  TinyVector<T,2> Data;
public:
  inline operator T() const     { return Data[ActiveCopy]; }
  //  inline operator T&()           { return Data[ActiveCopy]; }
  void operator= (const T &val) { Data[ActiveCopy] = val;  }
  inline void AcceptCopy()      { Data[OLDMODE] = Data[NEWMODE];       }
  inline void RejectCopy()      { Data[NEWMODE] = Data[OLDMODE];       }
  inline const T& operator[](int i) const
    { return Data[i]; }
  inline T& operator[](int i) 
    { return Data[i]; }


};


template <class T>
class Mirrored1DClass
{
private:
  TinyVector<Array<T,1>,2> Data;
public:
  inline void resize (int m)          {Data[OLDMODE].resize(m);Data[NEWMODE].resize(m);}
  inline int size() const             { return Data[NEWMODE].size();             }
  inline Array<T,1>& data() const     { return Data[ActiveCopy];           }
  inline Array<T,1>& data()           { return Data[ActiveCopy];           }
  inline operator Array<T,1>&() const { return Data[ActiveCopy];           }
  inline operator Array<T,1>&()       { return Data[ActiveCopy];           }
  inline T  operator()(int i) const   { return Data[ActiveCopy](i);        }
  inline T& operator()(int i)         { return Data[ActiveCopy](i);        }
  inline const Array<T,1>& operator[](int i) const 
  { return Data[i]; }
  inline Array<T,1>& operator[](int i) 
  { return Data[i]; }
  inline void AcceptCopy ()           { Data[OLDMODE] = Data[NEWMODE];     }
  inline void RejectCopy ()           { Data[NEWMODE] = Data[OLDMODE];     }
  inline void AcceptCopy (int i)
  { Data[OLDMODE](i) = Data[NEWMODE](i); }
  inline void RejectCopy (int i)
  { Data[NEWMODE](i) = Data[OLDMODE](i); }
  inline void AcceptCopy (int start, int end)
  { Data[OLDMODE](Range(start,end)) = Data[NEWMODE](Range(start,end)); }
  inline void RejectCopy (int start, int end)
  { Data[NEWMODE](Range(start,end)) = Data[OLDMODE](Range(start,end)); }
};


template <class T>
class Mirrored2DClass
{
  friend class PathClass;
protected:
  TinyVector<Array<T,2>,2> Data;
public:
  inline void resize (int m, int n)      
    { Data[OLDMODE].resize(m,n); Data[NEWMODE].resize(m,n); }
  inline void resizeAndPreserve(int m, int n)
    { Data[OLDMODE].resizeAndPreserve(m,n); Data[NEWMODE].resizeAndPreserve(m,n);}
  inline int rows() const                
    { return Data[NEWMODE].rows(); }
  inline int cols() const                
    { return Data[NEWMODE].cols(); }
  inline int extent(int i) const
    { return Data[NEWMODE].extent(i); }
  inline const Array<T,2>& data() const        
    { return Data[ActiveCopy]; }
  inline Array<T,2>& data()              
    { return Data[ActiveCopy]; }
  inline const Array<T,2>& operator[](int i) const
    { return Data[i]; }
  inline Array<T,2>& operator[](int i) 
    { return Data[i]; }
  inline operator Array<T,2>&() const    
    { return Data[ActiveCopy]; }
  inline operator Array<T,2>&()          
    { return Data[ActiveCopy]; }
  inline const T&  operator()(int i, int j) const 
    { return Data[ActiveCopy](i,j);      }
  inline T& operator()(int i, int j)       
    { return Data[ActiveCopy](i,j);      }
  inline void AcceptCopy ()              
    { Data[OLDMODE] = Data[NEWMODE]; }
  inline void RejectCopy ()              
    { Data[NEWMODE] = Data[OLDMODE]; }
};

template <class T>
class Mirrored3DClass
{
protected:
  TinyVector<Array<T,3>,2> Data;
public:
  inline void resize (int m, int n, int o)      
    { Data[NEWMODE].resize(m,n,o);  Data[OLDMODE].resize(m,n,o); }
  inline int rows() const                
    { return Data[NEWMODE].rows(); }
  inline int cols() const                
    { return Data[NEWMODE].cols(); }
  inline int depth() const 
    { return Data[NEWMODE].depth(); } 
  inline int extent (int i) const
    { return Data[NEWMODE].extent(i); }
  inline const Array<T,3>& data() const        
    { return Data[ActiveCopy]; }
  inline Array<T,3>& data()              
    { return Data[ActiveCopy]; }
  inline operator Array<T,3>&() const    
    { return Data[ActiveCopy]; }
  inline operator Array<T,3>&()          
    { return Data[ActiveCopy]; }
  inline const T& operator()(int i, int j, int k) const 
    { return Data[ActiveCopy](i,j,k);      }
  inline T& operator()(int i, int j, int k)       
    { return Data[ActiveCopy](i,j,k);      }
  inline const Array<T,3>& operator[](int i) const 
    { return Data[i]; }
  inline Array<T,3>& operator[](int i) 
    { return Data[i]; }
  inline void AcceptCopy ()              
    { Data[OLDMODE] = Data[NEWMODE]; }
  inline void RejectCopy ()              
    { Data[NEWMODE] = Data[OLDMODE]; }
};



// template <class T>
// class Mirrored2DClass
// {
// private:
//   TinyVector<Array<T,2>,2> Data;
// public:
//   inline void resize (int m, int n);
//   inline int size() const;
//   inline int rows() const;
//   inline int cols() const;
//   inline  T operator()(int i, i) const;
//   inline T& operator()(int i);
//   inline void AcceptCopy ();
//   inline void RejectCopy ();
//   inline void AcceptCopy (int start, int end);
//   inline void RejectCopy (int start, int end);
// };

void MirroredClassTest();
#endif
