/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "ReadPath.h"
///BUG: Doesn't load permutations yet
void 
ReadPathClass::MakeMove()
{
  //cerr << "Config " << filenames(fileIndex) << " " << currConfig << endl;
  Array<double,3> PathVec;
  IOVar->Read(PathVec,currConfig,Range::all(),Range::all(),Range::all());
  //cerr<<PathVec.extent(0)<<" "<<PathVec.extent(1)<<" "<<PathVec.extent(2)<<endl;
  //cerr<<PathData.Path.NumTimeSlices()<<endl;
  assert(PathVec.extent(0)==PathData.Path.NumParticles());
  if(PathVec.extent(1)!=PathData.Path.NumTimeSlices()-1) {
    if(overrideTimeSliceCheck) {
      //cerr << "Overriding time slice check; I will read the first " << PathData.Path.NumTimeSlices()-1 << " from a path with extent " << PathVec.extent(1) << endl;
    }
    else {
      cerr << "ERROR: Time slice mismatch " << PathData.Path.NumTimeSlices()-1 << " readpath has " << PathVec.extent(1) << "; ABORTING" << endl;
      exit(1);
    }
  }
  assert(PathVec.extent(2)==NDIM);
  
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
      dVec pos;
      for (int dim=0;dim<NDIM;dim++)
 	pos(dim)=PathVec(ptcl,slice,dim);
      PathData.Path.SetPos(slice,ptcl,pos);
      //   cerr<<PathVec.extent(0)<<" "<<PathVec.extent(1)<<" "<<PathVec.extent(2)<<endl;
    }
  }
  //cerr<<"Going to accept"<<endl;
  PathData.AcceptMove(0,PathData.Path.NumTimeSlices()-1,ActiveParticles);
  if(doAllConfigs) {
    currConfig++;
    if (currConfig>=NumConfigs){
      cerr<<"You are trying to read configurations that aren't there"<<endl;
      assert(1==2);
    }
  }
  else {
    configIndex++;
    if(configIndex == configIDs(fileIndex).size()) {
      cerr << "File " << fileIndex << " " << filenames(fileIndex) << " over.  Read " << configIndex << " configs." << endl;
      Close();
      configIndex = 0;
      fileIndex++;
      if(fileIndex == filenames.size()) {
        cerr << "All paths read.  Terminating." << endl;
        exit(1);
      }
      Init(filenames(fileIndex));
    }
    currConfig = configIDs(fileIndex)(configIndex);
  }
}

void ReadPathClass::Init(string filename)
{
  IOSectionClass in;
  assert (in.OpenFile(filename.c_str()));
  assert(in.OpenSection("Observables"));
  assert(in.OpenSection("PathDump"));
  IOVar = in.GetVarPtr("Path");
  Array<double,1> checkSize;
  IOVar->Read(checkSize,Range::all(),0,0,0);
  NumConfigs=checkSize.size();
  ActiveParticles.resize(PathData.Path.NumParticles());
  for (int i=0;i<ActiveParticles.size();i++)
    ActiveParticles(i)=i;
}

void ReadPathClass::Close() {
  delete IOVar;
  //IOVar.close();
  //in.CloseSection();
  //in.CloseSection();
  //in.Close();
}

void
ReadPathClass::Read(IOSectionClass &input)
{
  cerr<<"ReadPathClass::Read"<<endl;
  overrideTimeSliceCheck = false;
  input.ReadVar("OverrideTimeSlice",overrideTimeSliceCheck);
  string fileName;
  doAllConfigs = true;
  if(input.ReadVar("File",fileName)) {
    filenames.resize(1);
    filenames(0) = fileName;

  }
  else {
    cerr << "Using FileSection interface" << endl;
    doAllConfigs = false;
    assert(input.ReadVar("FileList",filenames));
    int numFiles = input.CountSections("FileSection");
    assert(numFiles == filenames.size());
    cerr << "Read " << numFiles << " files: " << filenames << endl;
    configIDs.resize(numFiles);
    for(int s=0; s<numFiles; s++) {
      input.OpenSection("FileSection",s);
      int numC;
      assert(input.ReadVar("NumConfigs",numC));
      configIDs(s).resize(numC);
      Array<int,1> cList;
      if(input.ReadVar("ConfigList",cList)) {
        assert(cList.size() == numC);
        for(int c=0; c<numC; c++)
          configIDs(s)(c) = cList(c);
      }
      else {
        int start, end, interval;
        assert(input.ReadVar("StartConfig",start));
        assert(input.ReadVar("EndConfig",end));
        assert(input.ReadVar("Interval",interval));
        int myC = start;
        int index = 0;
        cerr << "Got config parameters " << start << ", " << end << ", " << interval << endl;
        while(myC <= end) {
          configIDs(s)(index) = myC;
          myC += interval;
          index++;
        }
      }
      cerr << "Config list for file " << s << " is " << configIDs(s) << endl;
      input.CloseSection();
    }
  }

  fileIndex = 0;
  configIndex = 0;
  currConfig=0;
  Init(filenames(0));
}


