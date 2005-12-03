#ifndef OBSERVABLE_VAR_H
#define OBSERVABLE_VAR_H


#include "../Common.h"
#include "../PathDataClass.h"
#include <Common/IO/IO.h>

using namespace IO;

class ObservableVar
{
protected:
  IOSectionClass &Out;
  CommunicatorClass &Comm;
  bool FirstTime;
  IOVarBase *IOVar;
  string Name;
  
public:
  void Flush();
  ObservableVar (string name, IOSectionClass &out,  CommunicatorClass &comm) :
    Out(out), FirstTime(true), Comm(comm), Name(name)
  {
  }
};

class ObservableDouble : public ObservableVar
{
public:
  inline void Write (double val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<double,1> vec(1);
	vec(0) = val;
	Out.WriteVar (Name, vec);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableDouble (string name, IOSectionClass &out, CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};

class ObservableVecDouble1 : public ObservableVar
{
public:
  inline void Write (Array<double,1> &val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<double,2> mat(1,val.size());
	mat(0,Range::all()) = val;
	Out.WriteVar (Name, mat);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableVecDouble1(string name, IOSectionClass &out, 
		       CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};


 class ObservableVecInt1 : public ObservableVar
 {
 public:
   inline void Write (Array<int,1> &val)
   {
     if (Comm.MyProc()==0) {
       if (FirstTime) {
 	FirstTime=false;
 	Array<int,2> mat(1,val.size());
 	mat(0,Range::all()) = val;
 	Out.WriteVar (Name, mat);
 	IOVar = Out.GetVarPtr(Name);
       }
       else
 	IOVar->Append(val);
     }
   }
   ObservableVecInt1(string name, IOSectionClass &out, 
 		       CommunicatorClass &comm) 
     : ObservableVar (name, out, comm)
   {
     // do nothing
   }
 };

class ObservableVecDouble2 : public ObservableVar
{
public:
  inline void Write (Array<double,2> &val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<double,3> tensor(1,val.extent(0), val.extent(1));
	tensor(0,Range::all(),Range::all()) = val;
	Out.WriteVar (Name, tensor);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableVecDouble2(string name, IOSectionClass &out, 
		       CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};

class ObservableVecDouble3 : public ObservableVar
{
public:
  inline void Write (Array<double,3> &val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<double,4> tensor(1,val.extent(0), val.extent(1), val.extent(2));
	tensor(0,Range::all(),Range::all(),Range::all()) = val;
	Out.WriteVar (Name, tensor);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableVecDouble3(string name, IOSectionClass &out, 
		       CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};





class ObservableInt : public ObservableVar
{
public:
  inline void Write (int val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<int,1> vec(1);
	vec(0) = val;
	Out.WriteVar (Name, vec);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableInt (string name, IOSectionClass &out, CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};

class ObservableVecInt2 : public ObservableVar
{
public:
  inline void Write (Array<int,2> &val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<int,3> tensor(1,val.extent(0), val.extent(1));
	tensor(0,Range::all(),Range::all()) = val;
	Out.WriteVar (Name, tensor);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableVecInt2(string name, IOSectionClass &out, 
		       CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};

class ObservableVecInt3 : public ObservableVar
{
public:
  inline void Write (Array<int,3> &val)
  {
    if (Comm.MyProc()==0) {
      if (FirstTime) {
	FirstTime=false;
	Array<int,4> tensor(1,val.extent(0), val.extent(1), val.extent(2));
	tensor(0,Range::all(),Range::all(),Range::all()) = val;
	Out.WriteVar (Name, tensor);
	IOVar = Out.GetVarPtr(Name);
      }
      else
	IOVar->Append(val);
    }
  }
  ObservableVecInt3(string name, IOSectionClass &out, 
		       CommunicatorClass &comm) 
    : ObservableVar (name, out, comm)
  {
    // do nothing
  }
};






#endif
