#include <iostream>
#include <blitz/array.h> 
#include <hdf5.h>

using namespace std;

template<class T, int RANK>
class Var
{
protected:
  
public:
  virtual int rank() { return RANK; }
  virtual blitz::Array<T,RANK>& Read()                 = 0;
  virtual void Write(const blitz::Array<T,RANK> &data) = 0;
};

template<class T, int RANK>
class VarHDF5 : public Var<T, RANK>
{
private:
  typedef VarHDF5<T,RANK-1> SliceType;
  hid_t DataSetID, DataSpaceID, BoolType;
  const static hid_t AtomType;
public:
  int rank() { return 3; }
  blitz::Array<T,RANK>& Read();
  void Write(const blitz::Array<T,RANK> &data);
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
blitz::Array<T,RANK>& 
VarHDF5<T,RANK>::Read()
{

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


main()
{
  Var<double,4> *myVar;

  myVar = new VarHDF5<double,4>;

  cerr << "myVar.rank() = " << myVar->rank() << endl;
}
