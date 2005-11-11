#ifndef IO_VAR_BASE_H
#define IO_VAR_BASE_H

#include <blitz/array.h>
#include <hdf5.h>
#include <iostream>
using namespace blitz;

namespace IO {
  typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE, INVALID } IODataType;
  typedef enum { HDF5_TYPE, ASCII_TYPE} IOFileType;
  
  class IOVarBase
  {
  private:
    nilArraySection n0;
  protected:
    string Name;
  public:
    virtual int GetRank()            = 0;
    virtual IODataType GetType()     = 0;
    virtual IOFileType GetFileType() = 0;
    virtual void Print(std::ofstream &out) { /* do nothing */ }
    //    virtual string GetTypeString()   = 0;
    string GetName () const
    { return Name; }
    virtual int GetExtent(int dim)   = 0;

    /// Resizes the first dimension of the variable
    virtual void Resize(int n)       = 0;

    //////////////////////
    /// Read functions ///
    //////////////////////
    template<typename T> bool Read(T &val);
    template<typename T, int RANK> bool Read(Array<T,RANK> &val);
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9, typename T10>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8, T9 s9, T10 s10);
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8, T9 s9) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7) 
    { Rreturn ead(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { return Read(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { return Read(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { return Read(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { return Read(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { return Read(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1) 
    { return Read(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0>
    bool Read(Array<T,RANK> &val, T0 s0) 
    { return Read(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    

    ///////////////////////
    /// Write functions ///
    ///////////////////////    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9, typename T10>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8, T9 s9, T10 s10);
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8, typename T9>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8, T9 s9) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { return Write(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { return Write(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { return Write(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { return Write(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { return Write(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1) 
    { return Write(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0>
    bool Write(Array<T,RANK> &val, T0 s0) 
    { return Write(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
  
    template<typename T>
    bool Append (const T val)
    {

    }

    template<typename T, int RANK>
    bool Append(Array<T,RANK> &val)
    {
      assert (GetRank() == (RANK+1));
      int n = GetExtent(0);
      for (int i=0; i<RANK; i++)
	assert (val.extent(i) == GetExtent(i+1));
      Resize(n+1);
      return Write (val, n);
    }
    
  };


  ///////////////////////////////////////////////////////////////////
  /// The following are template tricks for converting a C++ type ///
  /// into an enumerated type variable.                           ///
  ///////////////////////////////////////////////////////////////////
  template<typename T> class TypeConvert
  { public: static const IODataType Type = INVALID; };

  template<> class TypeConvert<double>
  { public: static const IODataType Type = DOUBLE_TYPE; };

  template<> class TypeConvert<int>
  { public: static const IODataType Type = INT_TYPE; };

  template<> class TypeConvert<string>
  { public: static const IODataType Type = STRING_TYPE; };

  template<> class TypeConvert<bool>
  { public: static const IODataType Type = BOOL_TYPE; };
  



  template<typename T> inline string TypeString(T val) {  return "unknown"; }
  template<> string inline TypeString(double val)      {  return "double";  }
  template<> string inline TypeString(int val)         {  return "int";     }
  template<> string inline TypeString(string val)      {  return "string";  }
  template<> string inline TypeString(bool val)        {  return "bool";    }


}



#endif // ifndef IO_BASE_H
