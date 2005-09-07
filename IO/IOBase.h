#ifndef IO_BASE_H
#define IO_BASE_H

#include <blitz/array.h>
using namespace blitz;

typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE } IODataType;

class IOVarBase
{
public:
  virtual int GetRank()           = 0;
  virtual IODataType GetType()    = 0;
  virtual string GetTypeString()  = 0;
  template<typename T, int RANK> bool Read(Array<T,RANK> &val);
//   template<typename T,  typename T0, typename T1, typename T2,
// 	   typename T3, typename T4, typename T5, typename T6,
// 	   typename T7, typename T8, typename T9, typename T10>
//   bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val, 
// 	    T0 s0,T1 s1,T2 s2,T3 s3,T4 s4,T5 s5,T6 s6,T7 s7,T8 s8,T9
// s9,T10 s10);
  template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	   typename T3, typename T4, typename T5, typename T6, typename T7,
	   typename T8, typename T9, typename T10>
  bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	    T7 s7, T8 s8, T9 s9, T10 s10);
  
};


template<class T, int RANK>
class IOVar : public IOVarBase
{
protected:

public:
  int GetRank();
  DataType GetType();
  string GetTypeString();
  

};


// template<typename T,  typename T0, typename T1, typename T2, typename T3, 
// 	 typename T4, typename T5, typename T6, typename T7, typename T8, 
// 	 typename T9, typename T10> bool
// IOVarBase::Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
// 		T0 s0,T1 s1,T2 s2,T3 s3,T4 s4,T5 s5,T6 s6,T7 s7,T8 s8,T9 s9,T10 s10)
// {

// }	

template<typename T,  int RANK, typename T0, typename T1, typename T2, 
	 typename T3, typename T4, typename T5, typename T6, typename T7, 
	 typename T8, typename T9, typename T10> bool
IOVarBase::Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4,
		T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
{
    static const int numSlices = 
      SliceNum<T0>::isSlice+SliceNum<T1>::isSlice+SliceNum<T2>::isSlice+
      SliceNum<T3>::isSlice+SliceNum<T4>::isSlice+SliceNum<T5>::isSlice+
      SliceNum<T6>::isSlice+SliceNum<T7>::isSlice+SliceNum<T7>::isSlice+
      SliceNum<T8>::isSlice+SliceNum<T9>::isSlice+SliceNum<T10>::isSlice;
    
    /// The rank of the array must be the rank of the IO variable minus
    /// the number of slices by integer singlet ranges.
    static const int varRank=numSlices+RANK;

    if (((VarHDF5<T,varRank>*)newVar=
	 dynamic_cast<VarHDF5<T,varRank>*>(this)) != NULL) {
      
    }
    else if (((VarASCII<T,varRank>*)newVar=
	      dynamic_cast<VarASCII<T,varRank>*>(this)) != NULL) {

    }

}	

#endif // ifndef IO_BASE_H
