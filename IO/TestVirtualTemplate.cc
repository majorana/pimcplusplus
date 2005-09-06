#include <iostream>
#include <blitz/array.h> 
#include <hdf5.h>

using namespace std;
using namespace blitz;

typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE } DataType;

class VarBase
{
public:
  virtual int GetRank()       = 0;
  virtual DataType GetType()  = 0;
  virtual string GetTypeString() = 0;
};

template<class T, int RANK>
class Var : public VarBase
{
protected:
  
public:
  int GetRank();
  DataType GetType();
  string GetTypeString();
  virtual bool Read(blitz::Array<T,RANK>& val)                            = 0;
  virtual void Write(const blitz::Array<T,RANK> &data)                    = 0;

  virtual Var<T,RANK>& operator()(Range r0)                               = 0;
  virtual Var<T,RANK>& operator()(Range r0, Range r1)                     = 0;
  virtual Var<T,RANK>& operator()(Range r0, Range r1, Range r2)           = 0;
  virtual Var<T,RANK>& operator()(Range r0, Range r1, Range r2, Range r3) = 0;
  virtual Var<T,RANK>& operator()(Range r0, Range r1, Range r2, Range r3, 
				  Range r4)                               = 0;

};


template<class T, int RANK> int Var<T,RANK>::GetRank()
{  return RANK; }

///////////////////////
/// Double GetTypes ///
///////////////////////
template<> DataType Var<double,0>::GetType()
{  return DOUBLE_TYPE; }
template<> DataType Var<double,1>::GetType()
{  return DOUBLE_TYPE; }
template<> DataType Var<double,2>::GetType()
{  return DOUBLE_TYPE; }
template<> DataType Var<double,3>::GetType()
{  return DOUBLE_TYPE; }
template<> DataType Var<double,4>::GetType()
{  return DOUBLE_TYPE; }

template<> string Var<double,0>::GetTypeString()
{  return "double"; }
template<> string Var<double,1>::GetTypeString()
{  return "double"; }
template<> string Var<double,2>::GetTypeString()
{  return "double"; }
template<> string Var<double,3>::GetTypeString()
{  return "double"; }
template<> string Var<double,4>::GetTypeString()
{  return "double"; }


////////////////////
/// Int GetTypes ///
////////////////////
template<> DataType Var<int,0>::GetType()
{  return INT_TYPE; }
template<> DataType Var<int,1>::GetType()
{  return INT_TYPE; }
template<> DataType Var<int,2>::GetType()
{  return INT_TYPE; }
template<> DataType Var<int,3>::GetType()
{  return INT_TYPE; }
template<> DataType Var<int,4>::GetType()
{  return INT_TYPE; }

template<> string Var<int,0>::GetTypeString()
{  return "int"; }
template<> string Var<int,1>::GetTypeString()
{  return "int"; }
template<> string Var<int,2>::GetTypeString()
{  return "int"; }
template<> string Var<int,3>::GetTypeString()
{  return "int"; }
template<> string Var<int,4>::GetTypeString()
{  return "int"; }


///////////////////////
/// String GetTypes ///
///////////////////////
template<> DataType Var<string,0>::GetType()
{  return STRING_TYPE; }
template<> DataType Var<string,1>::GetType()
{  return STRING_TYPE; }
template<> DataType Var<string,2>::GetType()
{  return STRING_TYPE; }
template<> DataType Var<string,3>::GetType()
{  return STRING_TYPE; }
template<> DataType Var<string,4>::GetType()
{  return STRING_TYPE; }

template<> string Var<string,0>::GetTypeString()
{  return "string"; }
template<> string Var<string,1>::GetTypeString()
{  return "string"; }
template<> string Var<string,2>::GetTypeString()
{  return "string"; }
template<> string Var<string,3>::GetTypeString()
{  return "string"; }
template<> string Var<string,4>::GetTypeString()
{  return "string"; }



/////////////////////
/// Bool GetTypes ///
/////////////////////
template<> DataType Var<bool,0>::GetType()
{  return BOOL_TYPE; }
template<> DataType Var<bool,1>::GetType()
{  return BOOL_TYPE; }
template<> DataType Var<bool,2>::GetType()
{  return BOOL_TYPE; }
template<> DataType Var<bool,3>::GetType()
{  return BOOL_TYPE; }
template<> DataType Var<bool,4>::GetType()
{  return BOOL_TYPE; }

template<> string Var<bool,0>::GetTypeString()
{  return "bool"; }
template<> string Var<bool,1>::GetTypeString()
{  return "bool"; }
template<> string Var<bool,2>::GetTypeString()
{  return "bool"; }
template<> string Var<bool,3>::GetTypeString()
{  return "bool"; }
template<> string Var<bool,4>::GetTypeString()
{  return "bool"; }

class MyRange
{
public:
  int start, count, stride;
  int dim;
};


///////////////////////////////////////////////////////////////////////
///                       HDF5 Specialization                       ///
///////////////////////////////////////////////////////////////////////

template<typename T, int RANK> class VarHDF5;

template<typename T,  typename T0, typename T1, typename T2, typename T3, typename T4,  
         typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
class HDF5SliceMaker
{
public:
  static const int rank =      ArraySectionInfo<T0>::rank + ArraySectionInfo<T1>::rank + 
  ArraySectionInfo<T2>::rank + ArraySectionInfo<T3>::rank + ArraySectionInfo<T4>::rank + 
  ArraySectionInfo<T5>::rank + ArraySectionInfo<T6>::rank + ArraySectionInfo<T7>::rank + 
  ArraySectionInfo<T8>::rank + ArraySectionInfo<T9>::rank + ArraySectionInfo<T10>::rank;

  typedef VarHDF5<T,rank> SliceType;
};


template<class T, int RANK>
class VarHDF5 : public Var<T, RANK>
{
private:
  typedef VarHDF5<T,RANK> MyType;
  hid_t DatasetID, DiskSpaceID, MemSpaceID, BoolType;
  const static hid_t AtomType;
  TinyVector<MyRange,RANK> DiscRanges;
  TinyVector<int,RANK> Indices;

  template<typename T1, typename T2, typename T3, typename T4,  typename T5, typename T6,
	   typename T7, typename T8, typename T9, typename T10, typename T11>
  typename HDF5SliceMaker<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::SliceType 
  Slice(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6);

public:
  int rank() { return RANK; }
  
  /////////////////////////////////////////////////////////////////////
  /// The following Read functions read subarrays and slices of the ///
  /// current variable.                                             ///
  /////////////////////////////////////////////////////////////////////
  
  /// First, the real Read function which will be defined with real
  /// code with complete function specializations.
  bool Read (Array<T,RANK> &val);
  
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10>
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
	    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    return Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).Read(val);
  }
  
  template<typename T0> bool Read(typename SliceInfo<T,T0>::T_slice &val, 
				  T0 s0) 
  { nilArraySection n0;
    return Read(val, s0, n0, n0, n0, n0, n0, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1> bool Read(typename SliceInfo<T,T0,T1>::T_slice &val,
					       T0 s0, T1 s1) 
  { nilArraySection n0;
    return Read(val, s0, s1, n0, n0, n0, n0, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2> bool Read(typename SliceInfo<T,T0,T1,T2>::T_slice &val,
							    T0 s0, T1 s1, T2 s2) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, n0, n0, n0, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3>::T_slice &val, T0 s0, T1 s1, T2 s2, T3 s3) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, n0, n0, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4>::T_slice &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, n0, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5>::T_slice &val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, s5, n0, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6>::T_slice &val, 
	    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, s5, s6, n0, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7>::T_slice &val, 
	    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, n0, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8>::T_slice &val, 
	    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, n0, n0); }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9> 
  bool Read(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9>::T_slice &val, 
	    T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9) 
  { nilArraySection n0;
    return Read(val, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, n0); }

  ////////////////////////////////////////////////////////////////////////
  /// We use the Slice function to construct a subslice of the current ///
  /// variable.  This is called by the Read function.                  ///
  ////////////////////////////////////////////////////////////////////////
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10>
  typename HDF5SliceMaker<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType &
  Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10);

  void Write(const blitz::Array<T,RANK> &data);
  /// The following Write functions write subarrays and slices of the
  /// current variable.
  template<typename T1> 
  bool Write(typename SliceInfo<T,T1>::T_slice & val, T1 R1);
  template<typename T1, typename T2> 
  bool Write(typename SliceInfo<T,T1,T2>::T_slice & val, T1 R1, T2 R2);
  template<typename T1, typename T2, typename T3> 
  bool Write(typename SliceInfo<T,T1,T2,T3>::T_slice & val, T1 R1, T2 R2, T3 R3);
  template<typename T1, typename T2, typename T3, typename T4> 
  bool Write(typename SliceInfo<T,T1,T2,T3,T4>::T_slice & val, T1 R1, T2 R2, T3 R3, T4 R4);
  template<typename T1, typename T2, typename T3, typename T4, typename T5> 
  bool Write(typename SliceInfo<T,T1,T2,T3,T4,T5>::T_slice & val, T1 R1, T2 R2, T3 R3, T4 R4, T5 R5);

  /// Functions for writing a single atomic variable
  bool Write(T val, int i1);
  bool Write(T val, int i1, int i2);
  bool Write(T val, int i1, int i2, int i3);
  bool Write(T val, int i1, int i2, int i3, int i4);
  bool Write(T val, int i1, int i2, int i3, int i4, int i5);

  Var<T,RANK>& operator()(Range r0);
  Var<T,RANK>& operator()(Range r0, Range r1);
  Var<T,RANK>& operator()(Range r0, Range r1, Range r2);
  Var<T,RANK>& operator()(Range r0, Range r1, Range r2, Range r3);
  Var<T,RANK>& operator()(Range r0, Range r1, Range r2, Range r3, Range r4);
};

template<> const hid_t VarHDF5<double,0>::AtomType = H5T_NATIVE_DOUBLE;
template<> const hid_t VarHDF5<double,1>::AtomType = H5T_NATIVE_DOUBLE;
template<> const hid_t VarHDF5<double,2>::AtomType = H5T_NATIVE_DOUBLE;
template<> const hid_t VarHDF5<double,3>::AtomType = H5T_NATIVE_DOUBLE;
template<> const hid_t VarHDF5<double,4>::AtomType = H5T_NATIVE_DOUBLE;

template<> const hid_t VarHDF5<int,0>::AtomType = H5T_NATIVE_INT;
template<> const hid_t VarHDF5<int,1>::AtomType = H5T_NATIVE_INT;
template<> const hid_t VarHDF5<int,2>::AtomType = H5T_NATIVE_INT;
template<> const hid_t VarHDF5<int,3>::AtomType = H5T_NATIVE_INT;
template<> const hid_t VarHDF5<int,4>::AtomType = H5T_NATIVE_INT;


string TypeString (double a) { return ("double"); }
string TypeString (int    a) { return ("int");    }
string TypeString (string a) { return ("string"); }
string TypeString (bool   a) { return ("bool");   }
string TypeString (Range  a) { return ("Range");  }

template<class T, int RANK>
template<typename T1>
bool VarHDF5<T,RANK>::Write(typename SliceInfo<T,T1>::T_slice &val, T1 R1)
{
  T tval; T1 t1val;
  cerr << "Cannot index Array<" << TypeString(tval) << ", " << RANK
       << "> with a range of type " << TypeString(t1val) << endl;
  abort();
}

template<class T, int RANK>
bool VarHDF5<T,RANK>::Write(T val, int i1)
{
  T tval;
  cerr << "Cannot uniquely index a VarHDF5<" << TypeString(tval) << ", " << RANK 
       << "> with a single int..\n";
  abort();
}

template<class T, int RANK>
bool VarHDF5<T,RANK>::Write(T val, int i1, int i2)
{
  T tval;
  cerr << "Cannot uniquely index a VarHDF5<" << TypeString(tval) << ", " << RANK 
       << "> with a two integers.\n";
  abort();
}


template<class T, int RANK>
bool VarHDF5<T,RANK>::Write(T val, int i1, int i2, int i3)
{
  T tval;
  cerr << "Cannot uniquely index a VarHDF5<" << TypeString(tval) << ", " << RANK 
       << "> with a two integers.\n";
  abort();
}

template<>
bool VarHDF5<double,1>::Write(double val, int r1)
{
  cerr << "Writing range of HDF5 variable.\n";
}


template<>
template<>
bool VarHDF5<double,1>::Write<Range>(Array<double,1>& val, Range r1)
{
  cerr << "Writing range of HDF5 variable.\n";
}


template<>
template<>
bool VarHDF5<double,2>::Write<int,Range>(Array<double,1>& val, int r1, Range r2)
{
  cerr << "Writing range of HDF5 variable.\n";
}


template<>
template<>
bool VarHDF5<double,2>::Write<Range,int>(Array<double,1>& val, Range r1, int r2)
{
  cerr << "Writing range of HDF5 variable.\n";
}

template<>
template<>
bool VarHDF5<double,2>::Write<Range,Range>(Array<double,2>& val, Range r1, Range r2)
{
  cerr << "Writing range of HDF5 variable.\n";
}

////////////////////////////////////////////////////////////
/// Generic range operators:  These give error messages. ///
/// Specializations will allow valid cases.              ///
////////////////////////////////////////////////////////////
template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0)
{
  T a;
  cerr << "Var<" << TypeString(a) << "," << RANK 
       << "> cannot be ranged by (r0).\n";
  abort();
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1)
{
  T a;
  cerr << "Var<" << TypeString(a) << "," << RANK 
       << "> cannot be ranged by (r0, r1).\n";
  abort();

}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2)
{
  T a;
  cerr << "Var<" << TypeString(a) << "," << RANK 
       << "> cannot be ranged by (r0, r1).\n";
  abort();
}


template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2,
							Range r3)
{
  T a;
  cerr << "Var<" << TypeString(a) << "," << RANK 
       << "> cannot be ranged by (r0, r1).\n";
  abort();
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2,
							Range r3, Range r4)
{
  T a;
  cerr << "Var<" << TypeString(a) << "," << RANK 
       << "> cannot be ranged by (r0, r1).\n";
  abort();
}


/////////////////////////////////////////////////////////////////////
/// Specializations for range operator cases where they are valid ///
/////////////////////////////////////////////////////////////////////
template<>
Var<double,1>& VarHDF5<double,1>::operator()(Range r0)
{
  
}


template<class T, int RANK> bool
VarHDF5<T,RANK>::Read(blitz::Array<T,RANK>& val)
{
  cerr << "Reading variable of rank " << RANK << " and atomic size "
       << sizeof(T) << endl;
  return true;
}


template<> bool
VarHDF5<double,1>::Read(blitz::Array<double,1>& val)
{
  herr_t status = H5Dread(DatasetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &val);
  return (status == 0);
}


template<class T, int RANK>
void
VarHDF5<T,RANK>::Write(const blitz::Array<T,RANK> &data)
{

}

template<>
void
VarHDF5<double,4>::Write(const blitz::Array<double,4> &data)
{

}

template<class T, int RANK>
bool ReadVar(VarBase *var, blitz::Array<T,RANK> &val)
{
  string types[4] = {"double", "int", "string", "bool"};
  Var<T,RANK> *newVar = dynamic_cast<Var<T,RANK>*>(var);
  if (newVar == NULL) {
    T a;
    cerr << "Type mismatch in ReadVar:\n"
	 << "Variable has rank " << var->GetRank() 
	 << " and type " << var->GetTypeString() << ".\n"
	 << "Array    has rank " << RANK << " and type " << TypeString(a) 
	 << ".\n";
    abort();
  }
  else
    newVar->Read(val);
}


// template<class T, int RANK>
// bool VarHDF5<T,RANK>::Read(Array<T,RANK> &val, TinyVector<Range,RANK> ranges=Range::all())
// {
//   for (int i=0; i<RANK; i++) {
//     int index=Indices[i];
//     DataSetRanges(index) = ranges[i];
//   }
//   int n = DataSetRanges.size();
//   hsize_t start[n], count[n], stride[n];

//   for (int i=0; i<DataSetRanges.size(); i++) {
//     start(i) = DataSetRanges(i).first();
//     count(i) = DataSetRanges(i).length();
//     strid(i) = DataSetRanges(i).stride();
//   }
//   H5selectHyperslab(H5S_select, count, start, stride, NULL);
// }



template<class T, int RANK>
template<typename T1, typename T2, typename T3, typename T4,  typename T5, typename T6,
	 typename T7, typename T8, typename T9, typename T10, typename T11>
typename HDF5SliceMaker<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::SliceType 
VarHDF5<T,RANK>::Slice(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6)
{
  static const int numValid = SliceInfo<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::numValidTypes;
  assert(numValid <= RANK);
  int i=0;
  if(ArraySectionInfo<T1>::rank==0) {
  }
      
      

}

template<typename T>
class SubrangeType {
public:
  static const bool isInt=false, isRange=false, isNil=false;
};

template<>
class SubrangeType<Range> {
public:
  static const bool isInt=false, isRange=true, isNil=false;
};

template<>
class SubrangeType<int> {
public:
  static const bool isInt=true, isRange=false, isNil=false;
};

template<>
class SubrangeType<nilArraySection> {
public:
  static const bool isInt=false, isRange=false, isNil=false;
};


template<class T, int RANK>
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	 typename T6, typename T7, typename T8, typename T9, typename T10>
typename HDF5SliceMaker<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType &
VarHDF5<T,RANK>::Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
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
  
  newVar.MemSpaceID = H5Screate_simple(SliceInfo<T,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::rank,
				       memDims, memDims);

}





main()
{
  Var<double,4> *mydouble4Var;
  VarBase *myVar;
  
  mydouble4Var = new VarHDF5<double,4>;
  myVar = mydouble4Var;

  //  VarHDF5<double,3> &myVar2 = dynamic_cast<VarHDF5<double,3>&>
  //  (*myVar);
  Array<int,3> myArray;
  ReadVar(myVar, myArray);
  cerr << "myVar.GetRank() = " << myVar->GetRank() << endl;
  cerr << "myVar.GetType() = " << myVar->GetType() << endl;
}
