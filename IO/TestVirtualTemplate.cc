#include <iostream>
#include <blitz/array.h> 
#include <hdf5.h>

using namespace std;
using namespace blitz;

typedef enum { DOUBLE_TYPE, INT_TYPE, STRING_TYPE, BOOL_TYPE } DataType;

class VarBase
{
public:
  virtual int GetRank()      = 0;
  virtual DataType GetType() = 0;
};

template<class T, int RANK>
class Var : public VarBase
{
protected:
  
public:
  int GetRank();
  DataType GetType();
  virtual blitz::Array<T,RANK>& Read()                                    = 0;
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





template<class T, int RANK>
class VarHDF5 : public Var<T, RANK>
{
private:
  typedef VarHDF5<T,RANK> MyType;
  typedef VarHDF5<T,RANK-1> SliceType;
  hid_t DataSetID, DataSpaceID, BoolType;
  const static hid_t AtomType;
//   virtual void constructSubArray(VarHDF5<T,RANK>, Range r0);
//   virtual void constructSubArray(VarHDF5<T,RANK>, Range r0, Range r1);
//   virtual void constructSubArray(VarHDF5<T,RANK>, Range r0, Range r1,
// 				 Range r2);
//   virtual void constructSubArray(VarHDF5<T,RANK>, Range r0, Range r1,
// 				 Range r2, Range r3);

public:
  int rank() { return 3; }
  blitz::Array<T,RANK>& Read();
  void Write(const blitz::Array<T,RANK> &data);
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


template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0)
{
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1)
{
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2)
{
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2, 
					 Range r3)
{
}

template<class T, int RANK>
Var<T,RANK>& VarHDF5<T,RANK>::operator()(Range r0, Range r1, Range r2, 
					 Range r3, Range r4)
{
}


template<class T, int RANK>
blitz::Array<T,RANK>& 
VarHDF5<T,RANK>::Read()
{
  cerr << "Reading variable of rank " << RANK << " and atomic size "
       << sizeof(T) << endl;
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

string TypeString (double x)
{ return "double"; }
string TypeString (int x)
{ return "int"; }
string TypeString (string x)
{ return "string"; }
string TypeString (bool x)
{ return "bool"; }

template<class T, int RANK>
bool ReadVar(VarBase *var, blitz::Array<T,RANK> &val)
{
  string types[4] = {"double", "int", "string", "bool"};
  Var<T,RANK> *newVar = dynamic_cast<Var<T,RANK>*>(var);
  if (newVar == NULL) {
    T a;
    cerr << "Type mismatch in ReadVar:\n"
	 << "Variable has rank " << var->GetRank() 
	 << " and type " << types[var->GetType()] << ".\n"
	 << "Array    has rank " << RANK << " and type " << TypeString(a) 
	 << ".\n";
    abort();
  }
  else
    val = newVar->Read();
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
