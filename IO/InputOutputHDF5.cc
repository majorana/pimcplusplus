#include "InputOutputHDF5.h"
#include "InputOutputASCII.h"

/************************************************************
 *                    Input Functions                      *
 ************************************************************/

bool VarHDF5Class::ReadInto (double &val)
{
  assert (TypeClass == H5T_FLOAT);
  assert (Dimensions.size() == 1);
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &val);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<double,1> &val)
{
  assert (TypeClass == H5T_FLOAT);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<double,2> &val)
{
  assert (TypeClass == H5T_FLOAT);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<double,3> &val)
{
  assert (TypeClass == H5T_FLOAT);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}
 
bool VarHDF5Class::ReadInto (int &val)
{
  assert (TypeClass == H5T_INTEGER);
  assert (Dimensions.size() == 1);
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &val);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<int,1> &val)
{
  assert (TypeClass == H5T_INTEGER);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<int,2> &val)
{
  assert (TypeClass == H5T_INTEGER);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (Array<int,3> &val)
{
  assert (TypeClass == H5T_INTEGER);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (string &val)
{
  assert (TypeClass == H5T_STRING);
  assert (Dimensions.size() == 1);
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  Array<char,1> charArray(length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  val = charArray.data();
  H5Tclose(type);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (Array<string,1> &val)
{
  assert (TypeClass == H5T_STRING);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  Array<char,2> charArray(Dimensions(0),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    val(i) = &(charArray(i,0));
  H5Tclose(type);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (Array<string,2> &val)
{
  assert (TypeClass == H5T_STRING);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  Array<char,3> charArray(Dimensions(0),Dimensions(1),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    for (int j=0; j<Dimensions(1);j++)
      val(i,j) = &(charArray(i,j,0));
  H5Tclose(type);
  return (status == 0);
}



bool VarHDF5Class::ReadInto (Array<string,3> &val)
{
  assert (TypeClass == H5T_STRING);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  Array<char,4> charArray(Dimensions(0),Dimensions(1),
			  Dimensions(3),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    for (int j=0; j<Dimensions(1);j++)
      for (int k=0; k<Dimensions(2);k++)
	val(i,j,k) = &(charArray(i,j,k,0));
  H5Tclose(type);
  return (status == 0);
}



/// Strips everything after and including a '.' in the string.
/// Used to remove section numbers.
string InputTreeHDF5Class::StripName (string str)
{
  int pos = str.find(".");
  if (pos > 0)
    str.erase(pos,str.length());
  return(str);
}

bool InputTreeHDF5Class::OpenFile(string fileName,
				     string mySectionName,
				     InputTreeClass *parent)
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
  Name = mySectionName;
  return true;
}

/// C-style wrapper for member function iterator.
herr_t HDF5GroupIterator(hid_t group_id, const char *member_name,
			 void *classPtr)
{
  InputTreeHDF5Class &HDF5sec= *((InputTreeHDF5Class *)classPtr);
  HDF5sec.GroupIterator(member_name);
  return (0);
}

void InputTreeHDF5Class::GroupIterator(string member_name)
{
  //cerr << "GroupIterator( " << member_name << ")\n";

  H5G_stat_t statbuf;
  
  H5Gget_objinfo(GroupID, member_name.c_str(), 0, &statbuf);
  

  if (statbuf.type == H5G_GROUP) {
    InputTreeHDF5Class *newGroup = new InputTreeHDF5Class;
    newGroup->GroupID = H5Gopen (GroupID, member_name.c_str());
    if (newGroup->GroupID < 0) {
      cerr << "Error in InputTreeHDF5Class::GroupIterator.\n";
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
    int ndims = H5Sget_simple_extent_ndims(dataSpaceID);
    newVar->Dimensions.resize(ndims);
    H5Sget_simple_extent_dims(dataSpaceID, newVar->Dimensions.data(), 
			      NULL);
    H5Sclose (dataSpaceID);
    VarList.push_back(newVar);
  }
  else if (statbuf.type == H5G_TYPE) {
    cerr << "Compound types not yet supported "
	 << "in InputTreeHDF5Class.  Ignoring " 
	 << member_name << endl;
  }
  else
    cerr << " Unable to identify an object ";

}


void PrintIndent(int num)
{
  for (int counter=0;counter<num*3;counter++){
    cout<<' ';
  }
  

}



void InputTreeHDF5Class::PrintTree(int indentNum)
{
  PrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    PrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<InputTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

void InputTreeHDF5Class::PrintTree()
{
  PrintTree(0);
}

/// ReadGroup iterates over the members of it's group, creating
/// VarHDF5Class objects and new InputTreeHDF5Class objects as it
/// goes, calling itself recursively as necessary to traverse all the
/// subobjects below itself.
void InputTreeHDF5Class::ReadGroup(hid_t parentGroupID,
				      string name,
				      InputTreeClass *parent)
{
  Parent = parent;
  Name = StripName(name);  

  H5Giterate (parentGroupID, name.c_str(), (int *)NULL, HDF5GroupIterator,
	      this);
  // Make sure Iter is sane to start with.
  Iter = SectionList.begin();
}


void InputTreeHDF5Class::CloseFile()
{
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
   
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty())
    {
      SectionList.front()->CloseFile();
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



void OutputSectionHDF5Class::WriteVar(string name, int T)
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
				   H5T_NATIVE_INT, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}




void OutputSectionHDF5Class::WriteVar(string name, Array<int,1> &v)
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
				   H5T_NATIVE_INT, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}



void OutputSectionHDF5Class::WriteVar(string name, Array<int,2> &v)
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
				   H5T_NATIVE_INT, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}



void OutputSectionHDF5Class::WriteVar(string name, Array<int,3> &v)
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
				   H5T_NATIVE_INT, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   v.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}


void OutputSectionHDF5Class::WriteVar(string name,string str)
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
	    dim[i] = 1;
	  dataspace_id = H5Screate_simple(1, dim, NULL);
	  hid_t strType = H5Tcopy (H5T_C_S1);
	  H5Tset_size (strType, str.length()+1);
	    
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   strType, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, strType, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   str.c_str());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	  H5Tclose (strType);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}



void OutputSectionHDF5Class::WriteVar(string name,Array<string,1> &strs)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;

	  int maxLength = 0;
	  for (int i=0; i<strs.extent(0);i++)
	    if (strs(i).length() > maxLength)
		maxLength = strs(i).length();
	  maxLength++;
	  Array<char,2> charArray(strs.extent(0), maxLength);
	  for (int i=0; i<strs.extent(0); i++)
	    {
	      for (int x=0; x<strs(i).length(); x++)
		charArray(i,x) = (strs(i))[x];
	      int x = strs(i).length();
	      charArray(i,x) = '\0';
	    }

	  hsize_t dim[1];
	  for(int i=0; i<1; i++)
	    dim[i] = strs.extent(i);
	  dataspace_id = H5Screate_simple(1, dim, NULL);
	  hid_t strType = H5Tcopy (H5T_C_S1);
	  H5Tset_size (strType, maxLength);
	    
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   strType, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, strType, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   charArray.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	  H5Tclose (strType);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}





void OutputSectionHDF5Class::WriteVar(string name,Array<string,2> &strs)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;

	  int maxLength = 0;
	  for (int i=0; i<strs.extent(0);i++)
	    for (int j=0; j<strs.extent(1);j++)
	      if (strs(i,j).length() > maxLength)
		maxLength = strs(i,j).length();
	  maxLength++;
	  Array<char,3> charArray(strs.extent(0), strs.extent(1),
				  maxLength);
	  for (int i=0; i<strs.extent(0); i++)
	    for (int j=0; j<strs.extent(1);j++)
	      {
		for (int x=0; x<strs(i,j).length(); x++)
		  charArray(i,j,x) = (strs(i,j))[x];
		int x = strs(i,j).length();
		charArray(i,j,x) = '\0';
	    }

	  hsize_t dim[2];
	  for(int i=0; i<2; i++)
	    dim[i] = strs.extent(i);
	  dataspace_id = H5Screate_simple(2, dim, NULL);
	  hid_t strType = H5Tcopy (H5T_C_S1);
	  H5Tset_size (strType, maxLength);
	    
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   strType, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, strType, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   charArray.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	  H5Tclose (strType);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}




void OutputSectionHDF5Class::WriteVar(string name,Array<string,3> &strs)
{
  hid_t dataspace_id, dataset_id, loc_id;

  if (IsOpen)
    {
      if (SectionStack.empty()) 
	cerr << "Error in WriteVar:  No open sections.\n";
      else
	{
	  loc_id = SectionStack.top().GroupID;

	  int maxLength = 0;
	  for (int i=0; i<strs.extent(0);i++)
	    for (int j=0; j<strs.extent(1);j++)
	      for (int k=0; k<strs.extent(2);k++)
	      if (strs(i,j,k).length() > maxLength)
		maxLength = strs(i,j,k).length();
	  maxLength++;
	  Array<char,4> charArray(strs.extent(0), strs.extent(1),
				  strs.extent(2), maxLength);
	  for (int i=0; i<strs.extent(0); i++)
	    for (int j=0; j<strs.extent(1);j++)
	      for (int k=0; k<strs.extent(2);k++)
	      {
		for (int x=0; x<strs(i,j,k).length(); x++)
		  charArray(i,j,k,x) = (strs(i,j,k))[x];
		int x = strs(i,j,k).length();
		charArray(i,j,k,x) = '\0';
	    }

	  hsize_t dim[3];
	  for(int i=0; i<3; i++)
	    dim[i] = strs.extent(i);
	  dataspace_id = H5Screate_simple(3, dim, NULL);
	  hid_t strType = H5Tcopy (H5T_C_S1);
	  H5Tset_size (strType, maxLength);
	    
	  dataset_id =   H5Dcreate(loc_id, name.c_str(),
				   strType, dataspace_id,
				   H5P_DEFAULT);
	  herr_t status = H5Dwrite(dataset_id, strType, 
				   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
				   charArray.data());
	  if (status < 0)
	    cerr << "Error writing int to HDF5 file in WriteVar.\n";
	  H5Dclose (dataset_id);
	  H5Sclose (dataspace_id);
	  H5Tclose (strType);
	}
    }
  else
    cerr << "File not open in OutputSectionHDF5Class::WriteVar.\n";
}

