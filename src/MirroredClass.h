#ifndef MIRRORED_CLASS
#define MIRRORED_CLASS

#include <blitz/array.h>

extern int ActiveCopy;
typedef enum {OLDMODE, NEWMODE} ModeType;
inline void SetMode (ModeType type)
{
  if (type == OLDMODE)
    ActiveCopy = 0;
  else 
    ActiveCopy = 1;
}

using namespace blitz;

template <class T>
class MirroredClass
{
private:
  TinyVector<T,2> Data;
public:
  inline operator T() const     { return Data[ActiveCopy]; }
  void operator= (const T &val) { Data[ActiveCopy] = val;  }
  inline void AcceptCopy()      { Data[0] = Data[1];       }
  inline void RejectCopy()      { Data[1] = Data[0];       }
};


template <class T>
class Mirrored1DClass
{
private:
  TinyVector<Array<T,1>,2> Data;
public:
  inline void resize (int m)          {Data[0].resize(m);Data[1].resize(m);}
  inline int size() const             { return Data[0].size();             }
  inline Array<T,1>& data() const     { return Data[ActiveCopy];           }
  inline Array<T,1>& data()           { return Data[ActiveCopy];           }
  inline operator Array<T,1>&() const { return Data[ActiveCopy];           }
  inline operator Array<T,1>&()       { return Data[ActiveCopy];           }
  inline T  operator()(int i) const   { return Data[ActiveCopy](i);        }
  inline T& operator()(int i)         { return Data[ActiveCopy](i);        }
  inline void AcceptCopy ()           { Data[0] = Data[1];                 }
  inline void RejectCopy ()           { Data[1] = Data[0];                 }
  inline void AcceptCopy (int i)
  { Data[0](i) = Data[1](i); }
  inline void RejectCopy (int i)
  { Data[1](i) = Data[0](i); }
  inline void AcceptCopy (int start, int end)
  { Data[0](Range(start,end)) = Data[1](Range(start,end)); }
  inline void RejectCopy (int start, int end)
  { Data[1](Range(start,end)) = Data[0](Range(start,end)); }
};


template <class T>
class Mirrored2DClass
{
  friend class PathClass;
protected:
  TinyVector<Array<T,2>,2> Data;
public:
  inline void resize (int m, int n)      
    { Data[0].resize(m,n); Data[1].resize(m,n); }
  inline int rows() const                
    { return Data[0].rows(); }
  inline int cols() const                
    { return Data[0].cols(); }
  inline int extent(int i) const
    { return Data[0].extent(i); }
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
    { Data[0] = Data[1]; }
  inline void RejectCopy ()              
    { Data[1] = Data[0]; }
};

template <class T>
class Mirrored3DClass
{
protected:
  TinyVector<Array<T,3>,2> Data;
public:
  inline void resize (int m, int n, int o)      
    { Data[0].resize(m,n,o);  Data[1].resize(m,n,o); }
  inline int rows() const                
    { return Data[0].rows(); }
  inline int cols() const                
    { return Data[0].cols(); }
  inline int depth() const 
    { return Data[0].depth(); } 
  inline int extent (int i) const
    { return Data[0].extent(i); }
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
    { Data[0] = Data[1]; }
  inline void RejectCopy ()              
    { Data[1] = Data[0]; }
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
