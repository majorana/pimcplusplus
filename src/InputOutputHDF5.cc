#include "InputOutputHDF5.h"

bool OutputSectionHDF5Class::OpenFile(string fileName)
{
  FileName = fileName;
  bool success = true;
  
  FileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, 
		     H5P_DEFAULT, H5P_DEFAULT);
  success = IsOpen = (FileID >= 0);
  return (success);
}
      

void OutputSectionHD5Class::CloseFile()
{
  if (!GroupStack.empty())
    cerr << "You have not closed all open sections!  Error!\n";

  herr_t error = H5Fclose (FileID);
  if (error < 0)
    cerr << "Error closing file.\n";
}


void OutputSectionHDF5Class::OpenSection(string name)
{
  hid_t newGroup;
  if (GroupStack.empty())
    newGroup = H5Gcreate(FileID, name.c_str(), 0);
  else
    newGroup = H5Gcreate(GroupStack.top(), name.c_str(), 0);
  if (newGroup < 0)
    cerr << "Error adding group in OutputSectionHDF5Class::OpenSection"
	 << endl;
  else
    GroupStack.push (newGroup);
}

void OutputSectionHDF5Class::CloseSection()
{
  if (GroupStack.empty())
    cerr << "Error:  trying to close a non-existents section.\n";
  else
    {
      herr_t error = H5Gclose(GroupStack.top());
      if (error < 0)
	cerr << "Error closing section.\n";
      else
	GroupStack.pop();
    }
}



void OutputSectionHDF5Class::WriteVar(string name, double &T)
{
  hid_t dataspace_id, dataset_id, loc_id;


  if (GroupStack.empty()) // "root" variable
    loc_id = FileID;
  else
    loc_id = GroupStack.top();

  int dim[1] = 1;
  dataspace_id = H5Screate_simple(1, dim, NULL);
  dataset_id =   H5Dcreate(loc_id, name.c_str(),
			   H5T_INTEL_F64, dataspace_id,
			   H5P_DEFAULT);



}
