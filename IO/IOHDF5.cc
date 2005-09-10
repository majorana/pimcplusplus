#include "IOHDF5.h"

namespace IO {
  IOTreeClass *ReadTree (string fileName, string myName, IOTreeClass *parent);
  IOTreeClass *NewTree (string fileName, string myName, IOTreeClass *parent);

  IOFileType IOTreeHDF5Class::GetFileType()
  {
    return HDF5_TYPE;
  }

  /************************************************************
   *                     Helper Functions                     *
   ************************************************************/

  /// Strips everything after and including a '.' in the string.
  /// Used to remove section numbers.
  // void IOTreeHDF5Class::StripName (string str,string &newString,
  // 				    int &myInt)
  // {
  //   int pos = str.find(".");
  //   //  assert(pos>0);
  //   if (pos<=0){
  //     myInt=0;
  //     newString=str;
  //   }
  //   else{
  //     newString=str.substr(0,pos);
  //     string intString=str.substr(pos+1,str.length()-1);
  //     char *endptr;
  //     myInt=strtol(intString.c_str(),&endptr,10);
  //     assert (*endptr=='\0');
  //   }
  // }

  /// Strips everything after and including a '_' in the string.
  /// Used to remove section numbers.
  void IOTreeHDF5Class::StripName (string str, string &newString,
				   int &myInt)
  {
    int pos = str.length()-1;
    while ((pos>0) && (str[pos]!='_'))
      pos--;
    if (pos<=0){
      myInt=0;
      newString=str;
    }
    else{
      newString=str.substr(0,pos);
      string intString=str.substr(pos+1,str.length()-1);
      char *endptr;
      myInt=strtol(intString.c_str(),&endptr,10);
      assert (*endptr=='\0');
    }
  }


  /// Takes an integer and returns a string which is in the format .# 
  string NumExtension (int num)
  {
    string retString = "_";
    char numstr[100];
    snprintf(numstr, 100, "%d", num);
    //int len = strlen(numstr);
    //int numZeros = 8-len;
    //for (int i=0; i<numZeros; i++)
    //  retString += "0";
    retString += numstr;
    return (retString);
  }


  /************************************************************
   *                      File Functions                      *
   ************************************************************/

  bool IOTreeHDF5Class::OpenFile(string fileName,
				 string mySectionName,
				 IOTreeClass *parent)
  {
    cerr << "In HDF5 OpenFile.\n";
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto(&func, &client_data);
    H5Eset_auto(NULL, NULL);

    GroupID = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    BoolType = H5Topen(GroupID, "BOOL");
    // If we don't already have this in here (old version of this
    // library created this file), create a new one and stick it in the
    // file. 
    if (BoolType < 0) {
      BoolType = H5Tcreate(H5T_ENUM, sizeof(unsigned char));
      unsigned char val = 0;
      H5Tenum_insert(BoolType, "FALSE", &val);
      val = 1;
      H5Tenum_insert(BoolType, "TRUE", &val);
      H5Tcommit (GroupID, "BOOL", BoolType);
    }
    // And turn it back on;
    H5Eset_auto(func, client_data);


    if (GroupID < 0) {
      cerr << "Cannot open file " << fileName << endl;
      return false;
    }

    IsOpen = true;
    Parent = parent;
    ReadGroup (GroupID, "/", parent);
    Name = mySectionName;
    FileName=fileName;
    return true;
  }

  bool IOTreeHDF5Class::NewFile(string fileName,string myName,IOTreeClass* parent)
  {
    FileName = fileName;
    bool success = true;
  
    hid_t FileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, 
			     H5P_DEFAULT, H5P_DEFAULT);
    // Create a bool enum type
    BoolType = H5Tcreate(H5T_ENUM, sizeof(unsigned char));
    unsigned char val = (unsigned char) false;
    H5Tenum_insert(BoolType, "FALSE", &val);
    val = (unsigned char) true;
    H5Tenum_insert(BoolType, "TRUE", &val);
    H5Tcommit (FileID, "BOOL", BoolType);

    IsOpen = (FileID >= 0);
    success = IsOpen;
    if (success)
      {
	GroupID = FileID;
	Parent=parent;
	Name=myName;
      }
    return (success);
  }


  void IOTreeHDF5Class::CloseFile()
  {
    // First, free all the variables in the list
    while (!VarList.empty()) {
      delete(VarList.front());
      VarList.pop_front();
    }
   
    // Now, call all closes recursively and delete all sections
    while (!SectionList.empty()) {
      SectionList.front()->CloseFile();
      delete SectionList.front();
      SectionList.pop_front();
    }

    if (FileName!="") {
      // Release BoolType;
      H5Tclose (BoolType);
      H5Fclose(GroupID);
    }
    else {
      H5Gclose(GroupID);
    }
  }    


  void IOTreeHDF5Class::FlushFile()
  {
    // First, flush myself if I'm a root file.
    if (FileName!="") {
      herr_t err = H5Fflush(GroupID, H5F_SCOPE_GLOBAL);
      assert((int)err >= 0);
    }

    list<IOTreeClass*>::iterator iter = SectionList.begin();
    while (iter != SectionList.end()) {
      (*iter)->FlushFile();
      iter++;
    }
  }


  /************************************************************
   *                    Section Functions                     *
   ************************************************************/

  IOTreeClass* IOTreeHDF5Class::NewSection(string newName)
  {
    IOTreeHDF5Class* newSection = new IOTreeHDF5Class();
    newSection->Name=newName;
    newSection->Parent=this;
    newSection->MyNumber=CurrSecNum;
    newSection->BoolType = BoolType;
    SectionList.push_back(newSection);
  
    string numstr = NumExtension(CurrSecNum);
    newName += numstr;

    newSection->GroupID = H5Gcreate(GroupID,
				    newName.c_str(), 0);

    CurrSecNum++;
    return newSection;
  }

  void IOTreeHDF5Class::IncludeSection(IOTreeClass *newSection)
  {
    // Assign number
    newSection->MyNumber = CurrSecNum++;
    // Add to section list
    SectionList.push_back(newSection);
    // Put entry in my group in the HDF5 file
    string nameStr = newSection->Name + NumExtension (newSection->MyNumber);
    hid_t newGroupID = H5Gcreate(GroupID,
				 nameStr.c_str(), 0);
    hsize_t dim[1];
    dim[0] = 1;
    hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
    hid_t strType = H5Tcopy (H5T_C_S1);
    H5Tset_size (strType, newSection->FileName.length()+1);
    // Create a dataset named "include" which has the filename to include.
    hid_t dataset_id =   H5Dcreate(newGroupID, "include", strType, dataspace_id,
				   H5P_DEFAULT);
  
    herr_t status = H5Dwrite(dataset_id, strType, 
			     H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     newSection->FileName.c_str());
    if (status < 0)
      {
	cerr << "Error including to HDF5 file in IncludeSection.\n";
	exit(1);
      }

    H5Dclose (dataset_id);
    H5Sclose (dataspace_id);
    H5Tclose (strType);
    H5Gclose (newGroupID);
  }



  /// C-style wrapper for member function iterator.
  herr_t HDF5GroupIterator(hid_t group_id, const char *member_name,
			   void *classPtr)
  {
    IOTreeHDF5Class &HDF5sec= *((IOTreeHDF5Class *)classPtr);
    HDF5sec.GroupIterator(member_name);
    return (0);
  }


  void IOTreeHDF5Class::GroupIterator(string member_name)
  {
    //cerr << "GroupIterator( " << member_name << ")\n";

    H5G_stat_t statbuf;
  
    H5Gget_objinfo(GroupID, member_name.c_str(), 0, &statbuf);
  

    if (statbuf.type == H5G_GROUP) {
      // Open the group
      hid_t newGroupID;
      newGroupID = H5Gopen (GroupID, member_name.c_str());
      if ((int)newGroupID < 0) {
	cerr << "Error in IOTreeHDF5Class::GroupIterator.\n";
	exit(1);
      }
      // First check to see if we're including another file.
      // Turn off error printing
      H5E_auto_t func;
      void *client_data;
      H5Eget_auto(&func, &client_data);
      H5Eset_auto(NULL, NULL);
      hid_t includeData = H5Dopen(newGroupID, "include");
      // And turn it back on;
      H5Eset_auto(func, client_data);
      if ((int)includeData < 0) {
	// We're not including
	IOTreeHDF5Class *newTree = new IOTreeHDF5Class;
	newTree->GroupID = newGroupID;
	StripName(member_name,newTree->Name,newTree->MyNumber);
	InsertSection(newTree);
	newTree->BoolType = BoolType;
	newTree->ReadGroup (GroupID, member_name, this);
      }
      else {
	// We're including
	hid_t type = H5Dget_type(includeData);
	size_t length = H5Tget_size(type);
	blitz::Array<char, 1> charArray(length+1);
	herr_t status = H5Dread(includeData, type, H5S_ALL,
				H5S_ALL, H5P_DEFAULT, charArray.data());
	string fileName = charArray.data();
	H5Tclose(type);
	H5Dclose(includeData);
	string myName; int myNum;
	StripName (member_name, myName, myNum);
	IOTreeClass *newTree = ReadTree(fileName, myName, this);
	newTree->MyNumber = myNum;
	InsertSection (newTree);
      }
    }
    else if (statbuf.type == H5G_DATASET) {
      cerr << "Found HDF5 dataset:  " << member_name << endl;
      hid_t datasetID = H5Dopen(GroupID, member_name.c_str());
      VarList.push_back (NewIOVarHDF5(datasetID, member_name, BoolType));
    }
    else if (statbuf.type == H5G_TYPE) {
      if (member_name != "BOOL")
	cerr << "Compound types not yet supported "
	     << "in IOTreeHDF5Class.  Ignoring " 
	     << member_name << endl;
    }
    else
      cerr << " Unable to identify an object " << member_name << endl;;

  }


  /************************************************************
   *                    Printing Functions                    *
   ************************************************************/


  void PrintIndent(int num)
  {
    for (int counter=0;counter<num*2;counter++){
      cout<<' ';
    }
  }



  void IOTreeHDF5Class::PrintTree(int indentNum)
  {
    PrintIndent(indentNum);
    cout<<"Section: "<<Name<<endl;
    list<IOVarBase*>::iterator varIter=VarList.begin();
    while (varIter!=VarList.end()){
      PrintIndent(indentNum+1);
      cout<<"Variable: "<<(*varIter)->GetName()<<" "<<endl;
      varIter++;
    }
    list<IOTreeClass*>::iterator secIter=SectionList.begin();
    while (secIter!=SectionList.end()){
      //    cout<<"Section: "<<(*secIter)->Name<<endl;
      (*secIter)->PrintTree(indentNum+1);
      secIter++;
    }
  }

  void IOTreeHDF5Class::PrintTree()
  {
    PrintTree(0);
  }

  /// ReadGroup iterates over the members of it's group, creating
  /// VarHDF5Class objects and new IOTreeHDF5Class objects as it
  /// goes, calling itself recursively as necessary to traverse all the
  /// subobjects below itself.
  void IOTreeHDF5Class::ReadGroup(hid_t parentGroupID,
				  string name,
				  IOTreeClass *parent)
  {
    Parent = parent;
    int listLoc;
    StripName(name,Name,listLoc);  
  
    H5Giterate (parentGroupID, name.c_str(), (int *)NULL, HDF5GroupIterator,
		this);
    if (SectionList.empty())
      CurrSecNum = 0;
    else
      CurrSecNum = (SectionList.back())->MyNumber+1;
  }


}
