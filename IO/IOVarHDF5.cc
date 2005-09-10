#include "IOVarHDF5.h"

namespace IO {

  template<> bool
  IOVarHDF5<string,1>::VarRead(Array<string,1> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<string,2>::VarRead(Array<string,2> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<string,3>::VarRead(Array<string,3> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<string,4>::VarRead(Array<string,4> &val)
  {
    
  }
  
  
  template<> bool
  IOVarHDF5<bool,1>::VarRead(Array<bool,1> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<bool,2>::VarRead(Array<bool,2> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<bool,3>::VarRead(Array<bool,3> &val)
  {
    
  }
  
  template<> bool
  IOVarHDF5<bool,4>::VarRead(Array<bool,4> &val)
  {
    
  }
  
  
  IOVarBase *NewIOVarHDF5(hid_t dataSetID, string name, hid_t boolType)
  {
    /// First, figure out the rank
    hid_t diskSpaceID = H5Dget_space(dataSetID);
    hid_t memSpaceID  = H5Scopy(diskSpaceID);
    int rank = H5Sget_simple_extent_ndims(diskSpaceID);
    
    /// First, figure out what type it is.
    hid_t typeID = H5Dget_type(dataSetID);
    H5T_class_t classID = H5Tget_class(typeID);
    H5Tclose (typeID);
    if (classID == H5T_FLOAT) {
      if (rank == 1) 
	return new IOVarHDF5<double,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 2) 
	return new IOVarHDF5<double,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 3) 
	return new IOVarHDF5<double,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 4) 
	return new IOVarHDF5<double,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 5) 
	return new IOVarHDF5<double,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 6) 
	return new IOVarHDF5<double,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
    }
    if (classID == H5T_INTEGER) {
      if (rank == 1) 
	return  new IOVarHDF5<int,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 2) 
	return new IOVarHDF5<int,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 3) 
	return new IOVarHDF5<int,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 4) 
	return new IOVarHDF5<int,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 5) 
	return new IOVarHDF5<int,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 6) 
	return new IOVarHDF5<int,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
    }
    if (classID == H5T_STRING) {
      if (rank == 1) 
	return new IOVarHDF5<string,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 2) 
	return new IOVarHDF5<string,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 3) 
	return new IOVarHDF5<string,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 4) 
	return new IOVarHDF5<string,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 5) 
	return new IOVarHDF5<string,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 6) 
	return new IOVarHDF5<string,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
    }
    if (classID == H5T_ENUM){
      if (rank == 1) 
	return new IOVarHDF5<bool,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 2) 
	return new IOVarHDF5<bool,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 3) 
	return new IOVarHDF5<bool,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 4) 
	return new IOVarHDF5<bool,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 5) 
	return new IOVarHDF5<bool,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
      else if (rank == 6) 
	return new IOVarHDF5<bool,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, true);
    }
  }
}
