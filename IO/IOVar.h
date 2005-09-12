#ifndef IO_VAR_H
#define IO_VAR_H

#include "IOVarHDF5.h"
#include "IOVarASCII.h"

namespace IO {

  ///////////////////////////////////////////////////////////////////
  /// The following is a template trick for counting how many     ///
  /// dimension reductions we have made to a dataset by indexing  ///
  /// by integer arguements.                                      ///
  ///////////////////////////////////////////////////////////////////
  template<typename T> class SliceCheck
  { public:  static const int isSlice = 0; };

  template<> class SliceCheck<int>
  { public:  static const int isSlice = 1; };

  template<typename T, int RANK> bool
  IOVarBase::Read(Array<T,RANK> &val)
  {
    if (GetFileType() == HDF5_TYPE) {
      IOVarHDF5<T,RANK>* newVar = dynamic_cast<IOVarHDF5<T,RANK>*>(this);
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->VarRead(val);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,RANK>* newVar = dynamic_cast<IOVarASCII<T,RANK>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->VarRead(val);
    }
    else {
      cerr << "Error:  unknown type in IOVarBase::Read().\n";
      abort();
    }
  }	


  template<typename T,  int RANK, typename T0, typename T1, typename T2, 
	   typename T3, typename T4, typename T5, typename T6, typename T7, 
	   typename T8, typename T9, typename T10> bool
  IOVarBase::Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4,
		  T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    static const int numSlices = 
      SliceCheck<T0>::isSlice+SliceCheck<T1>::isSlice+SliceCheck<T2>::isSlice+
      SliceCheck<T3>::isSlice+SliceCheck<T4>::isSlice+SliceCheck<T5>::isSlice+
      SliceCheck<T6>::isSlice+SliceCheck<T7>::isSlice+SliceCheck<T7>::isSlice+
      SliceCheck<T8>::isSlice+SliceCheck<T9>::isSlice+SliceCheck<T10>::isSlice;
  
    /// The rank of the array must be the rank of the IO variable minus
    /// the number of slices by integer singlet ranges.
    static const int varRank=numSlices+RANK;
  
    if (GetFileType() == HDF5_TYPE) {
      IOVarHDF5<T,varRank>* newVar = dynamic_cast<IOVarHDF5<T,varRank>*>(this);
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->Slice(s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10).VarRead(val);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->Slice(s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10).VarRead(val);
    }
  }	


  template<typename T,  int RANK, typename T0, typename T1, typename T2, 
	   typename T3, typename T4, typename T5, typename T6, typename T7, 
	   typename T8, typename T9, typename T10> bool
  IOVarBase::Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4,
		   T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    static const int numSlices = 
      SliceCheck<T0>::isSlice+SliceCheck<T1>::isSlice+SliceCheck<T2>::isSlice+
      SliceCheck<T3>::isSlice+SliceCheck<T4>::isSlice+SliceCheck<T5>::isSlice+
      SliceCheck<T6>::isSlice+SliceCheck<T7>::isSlice+SliceCheck<T7>::isSlice+
      SliceCheck<T8>::isSlice+SliceCheck<T9>::isSlice+SliceCheck<T10>::isSlice;
    
    /// The rank of the array must be the rank of the IO variable minus
    /// the number of slices by integer singlet ranges.
    static const int varRank=numSlices+RANK;

    if (GetFileType() == HDF5_TYPE) {
      IOVarHDF5<T,varRank>* newVar = dynamic_cast<IOVarHDF5<T,varRank>*>(this);
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      return newVar->Slice(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10).VarWrite(val);
      //      return newVar->VarWriteSlice(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarASCII.\n";
	abort();
      }
      return newVar->Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarWrite(val);
    }

  }	

}       /// Ends namespace IO


#endif  /// Ends ifndef IO_VAR_H
