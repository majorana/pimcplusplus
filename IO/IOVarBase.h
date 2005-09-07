#ifndef IO_VAR_BASE_H
#define IO_VAR_BASE_H

#include <blitz/array.h>
#include <hdf5.h>
using namespace blitz;

namespace IO {
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
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7, T8 s8) 
    { Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	      T7 s7) 
    { Read(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { Read(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { Read(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { Read(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { Read(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { Read(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1>
    bool Read(Array<T,RANK> &val, T0 s0, T1 s1) 
    { Read(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0>
    bool Read(Array<T,RANK> &val, T0 s0) 
    { Read(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
    
    
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
    { Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7,
	     typename T8>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7, T8 s8) 
    { Write(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6, typename T7>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6,
	       T7 s7) 
    { Write(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5, typename T6>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
    { Write(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }
    
    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4, typename T5>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
    { Write(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3, typename T4>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
    { Write(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2,
	     typename T3>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2, T3 s3) 
    { Write(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1, typename T2>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1, T2 s2) 
    { Write(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0, typename T1>
    bool Write(Array<T,RANK> &val, T0 s0, T1 s1) 
    { Write(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }

    template<typename T,  int RANK,    typename T0>
    bool Write(Array<T,RANK> &val, T0 s0) 
    { Write(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }
  
  };


  ///////////////////////////////////////////////////////////
  ///                HDF5 Specializations                 ///  
  ///////////////////////////////////////////////////////////
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
    bool OwnDataset;
  public:
    int GetRank();
    IODataType GetType();
    IOFileType GetFileType();
    string GetTypeString();

    bool VarRead(T &val);
    bool VarRead(Array<T,RANK> &val);
    template<typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
		 T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);

    bool VarWrite(T &val);
    bool VarWrite(Array<T,RANK> &val);
    template<typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarWrite(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
		  T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);
    IOVarHDF5() : OwnDataset(false)
    {

    }
    
    IOVarHDF5(hid_t datasetID, hid_t diskSpaceID, hid_t memSpaceID, bool ownDataset=false)
    {
      DatasetID   = datasetID;
      DiskSpaceID = diskSpaceID;
      MemSpaceID  = memSpaceID;
      OwnDataset  = ownDataset;
    }
    ~IOVarHDF5()
    {
      H5Sclose(DiskSpaceID);
      H5Sclose(MemSpaceID);
      if (OwnDataset)
	H5Dclose(DatasetID);
      
    }

  };

  IOVarBase *NewHDF5Var(hid_t dataSetID)
  {
    /// First, figure out the rank
    hid_t diskSpaceID = H5Dget_space(dataSetID);
    hid_t memSpaceID  = H5Scopy(diskSpaceID);
    int rank = H5Sget_simple_extent_ndims(diskSpaceID);

    /// First, figure out what type it is.
    hid_t typeID = H5Dget_type(dataSetID);
    H5T_class_t classID = H5Tget_class(dataSetID);
    if (classID == H5T_FLOAT) {
      if (rank == 1) {
	IOVarHDF5<double,1> *newVar = new IOVarHDF5<double,1> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 2) {
	IOVarHDF5<double,2> *newVar = new IOVarHDF5<double,2> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 3) {
	IOVarHDF5<double,3> *newVar = new IOVarHDF5<double,3> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 4) {
	IOVarHDF5<double,4> *newVar = new IOVarHDF5<double,4> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 5) {
	IOVarHDF5<double,5> *newVar = new IOVarHDF5<double,5> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 6) {
	IOVarHDF5<double,6> *newVar = new IOVarHDF5<double,6> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
    }
    if (classID == H5T_INTEGER) {
      if (rank == 1) {
	IOVarHDF5<int,1> *newVar = new IOVarHDF5<int,1> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 2) {
	IOVarHDF5<int,2> *newVar = new IOVarHDF5<int,2> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 3) {
	IOVarHDF5<int,3> *newVar = new IOVarHDF5<int,3> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 4) {
	IOVarHDF5<int,4> *newVar = new IOVarHDF5<int,4> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 5) {
	IOVarHDF5<int,5> *newVar = new IOVarHDF5<int,5> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 6) {
	IOVarHDF5<int,6> *newVar = new IOVarHDF5<int,6> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
    }
    if (classID == H5T_STRING) {
      if (rank == 1) {
	IOVarHDF5<string,1> *newVar = new IOVarHDF5<string,1> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 2) {
	IOVarHDF5<string,2> *newVar = new IOVarHDF5<string,2> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 3) {
	IOVarHDF5<string,3> *newVar = new IOVarHDF5<string,3> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 4) {
	IOVarHDF5<string,4> *newVar = new IOVarHDF5<string,4> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 5) {
	IOVarHDF5<string,5> *newVar = new IOVarHDF5<string,5> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 6) {
	IOVarHDF5<string,6> *newVar = new IOVarHDF5<string,6> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
    }
    if (classID == H5T_ENUM){
      if (rank == 1) {
	IOVarHDF5<bool,1> *newVar = new IOVarHDF5<bool,1> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 2) {
	IOVarHDF5<bool,2> *newVar = new IOVarHDF5<bool,2> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 3) {
	IOVarHDF5<bool,3> *newVar = new IOVarHDF5<bool,3> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 4) {
	IOVarHDF5<bool,4> *newVar = new IOVarHDF5<bool,4> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 5) {
	IOVarHDF5<bool,5> *newVar = new IOVarHDF5<bool,5> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
      else if (rank == 6) {
	IOVarHDF5<bool,6> *newVar = new IOVarHDF5<bool,6> (dataSetID, diskSpaceID, memSpaceID, true);
	return newVar;
      }
    }


    H5Tclose(typeID);
  }

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

  // template<typename T,  typename T0, typename T1, typename T2, typename T3, 
  // 	 typename T4, typename T5, typename T6, typename T7, typename T8, 
  // 	 typename T9, typename T10> bool
  // IOVarBase::Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
  // 		T0 s0,T1 s1,T2 s2,T3 s3,T4 s4,T5 s5,T6 s6,T7 s7,T8 s8,T9 s9,T10 s10)
  // {

  // }	


  template<typename T, int RANK> 
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10> bool
  IOVarHDF5<T,RANK>::VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
			     T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
  {
    return Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarRead(val);
  }

  template<class T, int RANK>
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10>
  typename HDF5SliceMaker<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType &
  IOVarHDF5<T,RANK>::Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    typename HDF5SliceMaker<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType newVar;

    newVar.DatasetID = DatasetID;
    newVar.DiskSpaceID = H5Dget_space(DatasetID);
  
    hsize_t start[RANK], count[RANK], stride[RANK], dims[RANK], maxdims[RANK];
    hsize_t memDims[SliceInfo<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::rank];
    H5Sget_simple_extent_dims(newVar.DiskSpaceID, dims, maxdims);
  
    int memDimsIndex=0;
  
  
    /// Select the disk space hyperslab
    if (RANK > 0) {
      Range r0(s0);
      start[0]  = r0.first(0);
      count[0]  = (r0.last(dims[0]-1)-start[0])/r0.stride + 1;
      stride[0] = r0.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[0];
	memDimsIndex++;
      }
    }
    if (RANK > 1) {
      Range r1(s1);
      start[1] = r1.first(1);
      count[1] = (r1.last(dims[1]-1)-start[1])/r1.stride + 1;
      stride[1] = r1.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[1];
	memDimsIndex++;
      }
    }
    if (RANK > 2) {
      Range r2(s2);
      start[2] = r2.first(2);
      count[2] = (r2.last(dims[2]-1)-start[2])/r2.stride + 1;
      stride[2] = r2.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[2];
	memDimsIndex++;
      }
    }
    if (RANK > 3) {
      Range r3(s3);
      start[3] = r3.first(3);
      count[3] = (r3.last(dims[3]-1)-start[3])/r3.stride + 1;
      stride[3] = r3.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[3];
	memDimsIndex++;
      }
    }
    if (RANK > 4) {
      Range r4(s4);
      start[4] = r4.first(4);
      count[4] = (r4.last(dims[4]-1)-start[4])/r4.stride + 1;
      stride[4] = r4.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[4];
	memDimsIndex++;
      }
    }
    if (RANK > 5) {
      Range r5(s5);
      start[5] = r5.first(5);
      count[5] = (r5.last(dims[5]-1)-start[5])/r5.stride + 1;
      stride[5] = r5.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[5];
	memDimsIndex++;
      }
    }
    if (RANK > 6) {
      Range r6(s6);
      start[6] = r6.first(6);
      count[6] = (r6.last(dims[6]-1)-start[6])/r6.stride + 1;
      stride[6] = r6.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[6];
	memDimsIndex++;
      }
    }
    if (RANK > 7) {
      Range r7(s7);
      start[7] = r7.first(7);
      count[7] = (r7.last(dims[7]-1)-start[7])/r7.stride + 1;
      stride[7] = r7.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[7];
	memDimsIndex++;
      }
    }
    if (RANK > 8) {
      Range r8(s8);
      start[8] = r8.first(8);
      count[8] = (r8.last(dims[8]-1)-start[8])/r8.stride + 1;
      stride[8] = r8.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[8];
	memDimsIndex++;
      }
    }
    if (RANK > 9) {
      Range r9(s9);
      start[9] = r9.first(9);
      count[9] = (r9.last(dims[9]-1)-start[9])/r9.stride + 1;
      stride[9] = r9.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[9];
	memDimsIndex++;
      }
    }
    if (RANK > 10) {
      Range r10(s10);
      start[10] = r10.first(10);
      count[10] = (r10.last(dims[10]-1)-start[10])/r10.stride + 1;
      stride[10] = r10.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[10];
	memDimsIndex++;
      }
    }

    newVar.MemSpaceID = H5Screate_simple(SliceInfo<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::rank,
					 memDims, memDims);

  }


  template<typename T> class SliceCheck
  { public:  static const int isSlice = 0; };

  template<> class SliceCheck<int>
  { public:  static const int isSlice = 1; };


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
      newVar->VarRead(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      if (newVar == NULL) {
	cerr << "Error in dynamic cast to IOVarHDF5.\n";
	abort();
      }
      newVar->VarRead(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
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
      newVar->VarWrite(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
    }
    else if (GetFileType() == ASCII_TYPE) {
      IOVarASCII<T,varRank>* newVar = dynamic_cast<IOVarASCII<T,varRank>*>(this); 
      newVar->VarWrite(val, s0, s1, s2, s2, s4, s5, s6, s7, s8, s9, s10);
    }

  }	


  template<typename T> string TypeString(T val) {  return "unknown"; }
  template<> string TypeString(double val)      {  return "double";  }
  template<> string TypeString(int val)         {  return "int";     }
  template<> string TypeString(string val)      {  return "string";  }
  template<> string TypeString(bool val)        {  return "bool";    }

  /// This routine should cover double and int types.  Strings and bools
  /// need to be handled explicitly
  template<class T, int RANK> bool
  IOVarHDF5<T,RANK>::VarRead(Array<T,RANK> &val)
  {
    IODataType dataType = TypeConvert<T>::Type;
    hid_t memType;
    if      (dataType == DOUBLE_TYPE) memType = H5T_NATIVE_DOUBLE;
    else if (dataType == INT_TYPE)    memType = H5T_NATIVE_INT;
    else {
      T a;
      cerr << "Unknown data type in IOVarHDF5<" << TypeString(a) << ", " 
	   << RANK << ">" << endl;
    }
    /// Now, call HDF5 to do the actual reading.
    herr_t status = 
      H5Dread (DatasetID, memType, MemSpaceID, DiskSpaceID, H5P_DEFAULT, val.data());

  }

  template<> bool IOVarHDF5<string,0>::VarRead(string &val);
  template<> bool IOVarHDF5<string,1>::VarRead(Array<string,1> &val);
  template<> bool IOVarHDF5<string,2>::VarRead(Array<string,2> &val);
  template<> bool IOVarHDF5<string,3>::VarRead(Array<string,3> &val);
  template<> bool IOVarHDF5<string,4>::VarRead(Array<string,4> &val);
  template<> bool IOVarHDF5<bool,  0>::VarRead(bool &val);
  template<> bool IOVarHDF5<bool,  1>::VarRead(Array<bool,1> &val);
  template<> bool IOVarHDF5<bool,  2>::VarRead(Array<bool,2> &val);
  template<> bool IOVarHDF5<bool,  3>::VarRead(Array<bool,3> &val);
  template<> bool IOVarHDF5<bool,  4>::VarRead(Array<bool,4> &val);
  template<> bool IOVarHDF5<string,0>::VarWrite(string &val);
  template<> bool IOVarHDF5<string,1>::VarWrite(Array<string,1> &val);
  template<> bool IOVarHDF5<string,2>::VarWrite(Array<string,2> &val);
  template<> bool IOVarHDF5<string,3>::VarWrite(Array<string,3> &val);
  template<> bool IOVarHDF5<string,4>::VarWrite(Array<string,4> &val);
  template<> bool IOVarHDF5<bool,  0>::VarWrite(bool &val);
  template<> bool IOVarHDF5<bool,  1>::VarWrite(Array<bool,1> &val);
  template<> bool IOVarHDF5<bool,  2>::VarWrite(Array<bool,2> &val);
  template<> bool IOVarHDF5<bool,  3>::VarWrite(Array<bool,3> &val);
  template<> bool IOVarHDF5<bool,  4>::VarWrite(Array<bool,4> &val);
}



#endif // ifndef IO_BASE_H
