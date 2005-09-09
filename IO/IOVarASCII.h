#ifndef IO_VAR_ASCII_H
#define IO_VAR_ASCII_H

#include "IOVarBase.h"

namespace IO {

  template<typename T, int RANK>
  class IOVarASCII : public IOVarBase
  {
  protected:
    Array<T,RANK> MyValue;
  public:
    int GetRank();
    IODataType GetType();
    IOFileType GetFileType();
    string GetTypeString();

    int GetExtent(int dim);
    void Resize(int n);

    bool VarRead(Array<T,RANK> &val);
    template<typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
		 T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);

    bool VarWrite(Array<T,RANK> &val);
    template<typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarWrite(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
		  T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);
  };


  template<typename T, int RANK> int
  IOVarASCII<T,RANK>::GetRank()
  {
    return RANK;
  }

  template<typename T, int RANK> IODataType
  IOVarASCII<T,RANK>::GetType()
  {
    return TypeConvert<T>::IODataType;
  }

  template<typename T, int RANK> IOFileType
  IOVarASCII<T,RANK>::GetFileType()
  {
    return ASCII_TYPE;
  }


  template<typename T, int RANK> int
  IOVarASCII<T,RANK>::GetExtent(int dim) {
    return MyValue.extent(dim);
  }


  template<typename T, int RANK> void
  IOVarASCII<T,RANK>::Resize(int n) {
    TinyVector<int,RANK> dims = MyValue.shape();
    dims[0] = n;
    MyValue.resizeAndPreserve(dims);
  }

  template<typename T, int RANK> 
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10> bool
  IOVarASCII<T,RANK>::VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
			      T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
  {
    val = MyValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10);
  }


  template<typename T, int RANK> 
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10> bool
  IOVarASCII<T,RANK>::VarWrite(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
			       T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
  {
    MyValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10) = val;
  }

}


#endif
