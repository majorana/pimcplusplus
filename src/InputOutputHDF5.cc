#include "InputOutputHDF5.h"
#include "InputOutputASCII.h"

/************************************************************
 *                    Input Functions                      *
 ************************************************************/

bool InputSectionHDF5Class::OpenFile(string fileName, 
				     InputSectionClass *parent)
{
  GroupID = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (GroupID < 0)
    {
      cerr << "Cannot open file " << fileName << endl;
      return false;
    }

  IsOpen = true;
  IsRoot = true;
  Parent = parent;
  ReadGroup (GroupID, "/", NULL);
  Name = "all";
  return true;
}

/// C-style wrapper for member function iterator.
herr_t HDF5GroupIterator(hid_t group_id, const char *member_name,
			 void *classPtr)
{
  InputSectionHDF5Class &HDF5sec= *((InputSectionHDF5Class *)classPtr);
  HDF5sec.GroupIterator(member_name);
  return (0);
}

void InputSectionHDF5Class::GroupIterator(string member_name)
{
  cerr << "GroupIterator( " << member_name << ")\n";

  H5G_stat_t statbuf;
  
  H5Gget_objinfo(GroupID, member_name.c_str(), 0, &statbuf);
  

  if (statbuf.type == H5G_GROUP) {
    InputSectionHDF5Class *newGroup = new InputSectionHDF5Class;
    newGroup->GroupID = H5Gopen (GroupID, member_name.c_str());
    if (newGroup->GroupID < 0) {
      cerr << "Error in InputSectionHDF5Class::GroupIterator.\n";
      exit(10);
    }
    SectionList.push_back(newGroup);
    newGroup->ReadGroup (GroupID, member_name, this);
  }
  else if (statbuf.type == H5G_DATASET) {
    VarHDF5Class *newVar = new VarHDF5Class;
    newVar->DataSetID = H5Dopen(GroupID, member_name.c_str());
    newVar->Name = member_name;
    hid_t dataTypeID = H5Dget_type (newVar->DataSetID);
    newVar->TypeClass = H5Tget_class (dataTypeID);
    H5Tclose (dataTypeID);
    hid_t dataSpaceID = H5Dget_space (newVar->DataSetID);
    newVar->Ndims = H5Sget_simple_extent_ndims(dataSpaceID);
    newVar->Dimensions.resize(newVar->Ndims);
    H5Sget_simple_extent_dims(dataSpaceID, newVar->Dimensions.data(), NULL);
    H5Sclose (dataSpaceID);
    VarList.push_back(newVar);
  }
  else if (statbuf.type == H5G_TYPE) {
    cerr << "Compound types not yet supported "
	 << "in InputSectionHDF5Class.  Ignoring " 
	 << member_name << endl;
  }
  else
    cerr << " Unable to identify an object ";

}

/// ReadGroup iterates over the members of it's group, creating
/// VarHDF5Class objects and new InputSectionHDF5Class objects as it
/// goes, calling itself recursively as necessary to traverse all the
/// subobjects below itself.
void InputSectionHDF5Class::ReadGroup(hid_t parentGroupID,
				      string name,
				      InputSectionClass *parent)
{
  Parent = parent;
  Name = name;
  
  cerr << "name = " << name << endl;

  H5Giterate (parentGroupID, name.c_str(), (int *)NULL, HDF5GroupIterator,
	      this);


}


void InputSectionHDF5Class::Close()
{
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
   
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty())
    {
      SectionList.front()->Close();
      delete SectionList.front();
      SectionList.pop_front();
    }
  if (IsRoot)
    H5Fclose(GroupID);
  else
    H5Gclose(GroupID);
}    



/************************************************************
 *                    Output Functions                      *
 ************************************************************/ 

bool OutputSectionHDF5Class::OpenFile(string fileName)
{
  FileName = fileName;
  bool success = true;
  
  FileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, 
		     H5P_DEFAULT, H5P_DEFAULT);
  IsOpen = (FileID >= 0);
  success = IsOpen;
  if (success)
    {
      HDF5SectionClass newSection;
      newSection.GroupID = FileID;
      SectionStack.push(newSection);
    }
  return (success);
}
      

void OutputSectionHDF5Class::CloseFile()
{
  SectionStack.pop();
  if (!SectionStack.empty())
    cerr << "You have not closed all open sections!  Error!\n";

  herr_t error = H5Fclose (FileID);
  if (error < 0)
    cerr << "Error closing file.\n";
}


void OutputSectionHDF5Class::OpenSection(string name)
{
  hid_t newGroup;
  int num;
  HDF5SectionClass newSection;  

  if (SectionStack.empty())
    {
      num = newSection.SectionNumber(name);
      char numstr[100];
      snprintf (numstr, 100, ".%d", num);
      name += numstr;
      newSection.GroupID = H5Gcreate(FileID, name.c_str(), 0);
    }
  else
    {
      num = SectionStack.top().SectionNumber(name);
      char numstr[100];
      snprintf (numstr, 100, ".%d", num);
      name += numstr;

      newSection.GroupID = H5Gcreate(SectionStack.top().GroupID,
				     name.c_str(), 0);
    }
  if (newSection.GroupID < 0)
    cerr << "Error adding group in OutputSectionHDF5Class::OpenSection"
	 << endl;
  else
      SectionStack.push (newSection);
}

void OutputSectionHDF5Class::CloseSection()
{
  if (SectionStack.empty())
    cerr << "Error:  trying to close a non-existents section.\n";
  else
    {
      herr_t error = H5Gclose(SectionStack.top().GroupID);
      if (error < 0)
	cerr << "Error closing section.\n";
      else
	SectionStack.pop();
    }
}



void OutputSectionHDF5Class::WriteVar(string name, double T)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;
	  
	  hsize_t dim[1];
	  dim[0] = 1;
	  dataspace_id = H5Screate_simple(1, dim, NULL);
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
	  if (status < 0)
	    cerr << "Error writing double to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}




void OutputSectionHDF5Class::WriteVar(string name, Array<double,1> &v)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;
	  
	  hsize_t dim[1];
	  for(int i=0; i<1; i++)
	    dim[i] = v.extent(i);
	  dataspace_id = H5Screate_simple(1, dim, NULL);
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing double to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}



void OutputSectionHDF5Class::WriteVar(string name, Array<double,2> &v)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;
	  
	  hsize_t dim[2];
	  for(int i=0; i<2; i++)
	    dim[i] = v.extent(i);
	  dataspace_id = H5Screate_simple(2, dim, NULL);
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing double to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}



void OutputSectionHDF5Class::WriteVar(string name, Array<double,3> &v)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;
	  
	  hsize_t dim[3];
	  for(int i=0; i<3; i++)
	    dim[i] = v.extent(i);
	  dataspace_id = H5Screate_simple(3, dim, NULL);
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing double to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}


