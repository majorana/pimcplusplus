#ifndef IO_BASE_H
#define IO_BASE_H

#include <blitz/array.h>
#include <hdf5.h>
using namespace blitz;

typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE, INVALID } IODataType;
typedef enum { HDF5_TYPE, ASCII_TYPE} IOFileType;

class IOVarBase
{
private:
  nilArraySection n0;
public:
  virtual int GetRank()            = 0;
  virtual IODataType GetType()     = 0;
  virtual IOFileType GetFileType() = 0;
  virtual string GetTypeString()   = 0;
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
  template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	   typename T3, typename T4, typename T5, typename T6, typename T7,
	   typename T8, typename T9>
  bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	    T7 s7, T8 s8, T9 s9) 
  { Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
  
};


template<typename T, int RANK> class IOVarHDF5;

template<typename T,  typename T0, typename T1, typename T2, typename T3, typename T4,  
         typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
class HDF5SliceMaker
{
public:
  static const int rank =      ArraySectionInfo<T0>::rank + ArraySectionInfo<T1>::rank + 
  ArraySectionInfo<T2>::rank + ArraySectionInfo<T3>::rank + ArraySectionInfo<T4>::rank + 
  ArraySectionInfo<T5>::rank + ArraySectionInfo<T6>::rank + ArraySectionInfo<T7>::rank + 
  ArraySectionInfo<T8>::rank + ArraySectionInfo<T9>::rank + ArraySectionInfo<T10>::rank;

  typedef IOVarHDF5<T,rank> SliceType;
};

template<typename T, int RANK>
class IOVarHDF5 : public IOVarBase
{
protected:
  hid_t DatasetID, DiskSpaceID, MemSpaceID;
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10>
  typename HDF5SliceMaker<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType &
  Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10);
public:
  int GetRank();
  IODataType GetType();
  IOFileType GetFileType();
  string GetTypeString();
  template<typename T0, typename T1, typename T2, typename T3, typename T4,
	   typename T5, typename T6, typename T7, typename T8, typename T9,
	   typename T10>
  bool VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
	       T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);

};


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
  template<typename T0, typename T1, typename T2, typename T3, typename T4,
	   typename T5, typename T6, typename T7, typename T8, typename T9,
	   typename T10>
  bool VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
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




template<typename T, int RANK> 
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	 typename T6, typename T7, typename T8, typename T9, typename T10> bool
IOVarASCII<T,RANK>::VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
			    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
{
  val = MyValue(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10);
}

// template<typename T,  typename T0, typename T1, typename T2, typename T3, 
// 	 typename T4, typename T5, typename T6, typename T7, typename T8, 
// 	 typename T9, typename T10> bool
// IOVarBase::Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
// 		T0 s0,T1 s1,T2 s2,T3 s3,T4 s4,T5 s5,T6 s6,T7 s7,T8 s8,T9 s9,T10 s10)
// {

// }	




template<typename T>
class SliceCheck
{
public:
  static const int isSlice = 0;
};

template<>
class SliceCheck<int>
{
public:
  static const int isSlice = 1;
};


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
      
      newVar->VarRead(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      newVar->VarRead(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
    }

}	

#endif // ifndef IO_BASE_H
