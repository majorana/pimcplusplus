#include "InputOutputHDF5.h"

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
      SectionClass newSection;
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
  SectionClass newSection;  

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
      if (SectionStack.empty()) // "root" variable
	loc_id = FileID;
      else
	loc_id = SectionStack.top().GroupID;
      
      hsize_t dim[1];
      dim[0] = 1;
      dataspace_id = H5Screate_simple(1, dim, NULL);
      dataset_id =   H5Dcreate(loc_id, name.c_str(),
			       H5T_INTEL_F64, dataspace_id,
			       H5P_DEFAULT);
      herr_t status = H5Dwrite(dataset_id, H5T_INTEL_F64, 
			       H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
      if (status < 0)
	cerr << "Error writing double to HDF5 file in WriteVar.\n";
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}
