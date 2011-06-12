#include <fstream>
#include <string.h>
#include <iostream>
#include "Blitz.h"
#include <list>

using namespace std;




class VarClass 
{  
public:
  string Name;

  virtual bool ReadInto (double &val) = 0;
  /*  virtual bool ReadInto (Array<double,1> &val) = 0;
  virtual bool ReadInto (Array<double,2> &val)     = 0;
  virtual bool ReadInto (Array<double,3> &val)     = 0;
  virtual bool ReadInto (int &val)                 = 0;
  virtual bool ReadInto (Array<int,1> &val)        = 0;
  virtual bool ReadInto (Array<int,2> &val)        = 0;
  virtual bool ReadInto (Array<int,3> &val)        = 0;
  virtual bool ReadInto (bool &val)                = 0;
  virtual bool ReadInto (Array<bool,1> &val)       = 0;
  virtual bool ReadInto (Array<bool,2> &val)       = 0;
  virtual bool ReadInto (Array<bool,3> &val)       = 0;
  virtual bool ReadInto (string &val)              = 0;
  virtual bool ReadInto (Array<string,1> &val)     = 0;
  virtual bool ReadInto (Array<string,2> &val)     = 0;
  virtual bool ReadInto (Array<string,3> &val)     = 0; */
};



class VarASCIIClass : public VarClass
{  
public:
  Array<char,1> Value;
  bool ReadInto (double &val)
  {
    int l = Value.size();
    char str[l+1];
    for (int i=0; i<l; i++)
      str[i] = Value(i);
    str[l] ='\0';
    val = atof(str);
    return true;
  }
  /*  bool ReadInto (Array<double,1> &val);
  bool ReadInto (Array<double,2> &val);
  bool ReadInto (Array<double,3> &val);
  bool ReadInto (int &val);
  bool ReadInto (Array<int,1> &val);
  bool ReadInto (Array<int,2> &val);
  bool ReadInto (Array<int,3> &val);
  bool ReadInto (bool &val);
  bool ReadInto (Array<bool,1> &val);
  bool ReadInto (Array<bool,2> &val);
  bool ReadInto (Array<bool,3> &val);
  bool ReadInto (string &val);
  bool ReadInto (Array<string,1> &val);
  bool ReadInto (Array<string,2> &val);
  bool ReadInto (Array<string,3> &val); */
};
inline bool checkPair(Array<char,1> &buffer,int counter,char* toSee)
{
  if (counter+1>=buffer.size()){
    return false;
  }
  if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
    return true;
  }
  else return false;

}


void ReadWithoutComments(string fileName,Array<char,1> &buffer)
{
  ifstream infile;
  infile.open(fileName.c_str());
  Array<char,1> tmpBuffer;
  int counter=0;
  bool inQuote=false;
  char dummyChar;
  while (!infile.eof()){
    infile.get(dummyChar);    
    counter++;
  }

  tmpBuffer.resize(counter);
  buffer.resize(counter);
  counter=-1;
  infile.close();
  ifstream infile2;
  infile2.open(fileName.c_str());
  while (!infile2.eof()){
    counter++;
    infile2.get(dummyChar);
    tmpBuffer(counter)=dummyChar;    
  }
  //  cout<<tmpBuffer;
  int bufferLoc=0;
  for (int counter=0;counter<tmpBuffer.size();counter++){
    if (inQuote){
      while (tmpBuffer(counter) != '\"' && counter<tmpBuffer.size()){
	buffer(bufferLoc)=tmpBuffer(counter);
	counter++;
	bufferLoc++;
      }
      buffer(bufferLoc)=tmpBuffer(counter);
      bufferLoc++;
      inQuote=false;
    }
    else {
      if (checkPair(tmpBuffer,counter,"//")){
	while (tmpBuffer(counter)!='\n' && counter<tmpBuffer.size()){
	  counter++;
	}
	buffer(bufferLoc)=tmpBuffer(counter); //copy the \n over
	bufferLoc++;
      }
      else if (checkPair(tmpBuffer,counter,"/*")){
	while (!checkPair(tmpBuffer,counter,"*/") && counter<tmpBuffer.size()){
	  counter++;
	}
	counter++; //end up in the / of comment
      }
      else if (tmpBuffer(counter)=='\"'){
	inQuote=true;
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
      else {
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
    }
  }
  buffer.resizeAndPreserve(bufferLoc+1);
}



class SectionClass
{
public:
  SectionClass* parent;
  list<SectionClass*> SectionList;
  list<VarClass*> VarList;
  list<SectionClass*>::iterator Iter;
public:
  string Name;
  bool FindSection (string name, SectionClass * &sectionPtr, 
		    bool rewind=true);
  template<class T> bool ReadVar(string name, T &var);
  virtual bool ReadFile (string fileName) = 0;

};

bool SectionClass::FindSection (string name, SectionClass* &sectionPtr,
				    bool rewind)
{

  while ((Iter != SectionList.end()) && ((*Iter)->Name!=name)){
    Iter++;
  }
  bool found = Iter != SectionList.end(); 
  if (found)
    sectionPtr = *Iter;
  if (rewind)
    Iter = SectionList.begin();
  return (found);
}


template <class vartype>
class SpecialSectionClass : public SectionClass
{
public:
  bool ReadFile (string fileName);
};



template<class T>
bool SectionClass::ReadVar(string name, T &var)
{
  bool readVarSuccess;
  list<VarClass*>::iterator varIter=VarList.begin();
  while ((varIter!=VarList.end() && (*varIter)->Name!=name)){
    varIter++;
  }
  bool found = varIter != VarList.end();
  if (found)
    readVarSuccess=(*varIter)->ReadInto(var);
  else if (parent!=NULL){
    readVarSuccess=parent->ReadVar(name,var);
  }
  else {
    cerr<<"Couldn't find variable "<<name<<endl;
    cerr<<"Why am I getting called!"<<endl;
    return false;
  }  
  return readVarSuccess;	 
}



// Dummy:

bool isNumChar(char theChar)
{
  if ((theChar>='A' && theChar<='Z') ||
      (theChar>='a' && theChar<='z') ||
      (theChar>='1' && theChar<='9') ||
      theChar=='_' ||
      theChar=='\"'||
      theChar=='-')
    return true;
  return false;
}
void printCharArray(Array <char,1> theArray)
{
  for (int counter=0;counter<theArray.size();counter++){
    cout<<theArray(counter);
  }
  cout<<endl;
}

void getVariable(Array <char,1> &buffer,int &counter,
		 string &theName,Array <char,1> &theValue)
{
  theName="";
  int eqLocation=counter;
  while ((!(isNumChar(buffer(counter)))) &&
	 counter>0){
    counter--;
  }

  int endLocInclusive=counter;
  while ((isNumChar(buffer(counter)))
	 && counter>0){
    counter--;
  }
  int beginLocInclusive=counter+1;
  
  for (int strCount=beginLocInclusive;strCount<=endLocInclusive;
       strCount++){
    theName=theName+buffer(strCount);
  }
  ///Got name..now get value
  bool inQuotes=false;
  beginLocInclusive=eqLocation+1;
  while (!isNumChar(buffer(beginLocInclusive))){
    beginLocInclusive++;
  }
  endLocInclusive=beginLocInclusive;
  while (isNumChar(buffer(endLocInclusive)) ||  inQuotes){
    if (buffer(endLocInclusive)=='\"' && buffer(endLocInclusive-1)!='\\')
      inQuotes=!inQuotes;    
    endLocInclusive++;
  }
  endLocInclusive=endLocInclusive-1;
  theValue.resize(endLocInclusive-beginLocInclusive+1);
  int strCount2=0;
  for (int strCount=beginLocInclusive;strCount<=endLocInclusive;
       strCount++){
    theValue(strCount2)=buffer(strCount);
    strCount2++;
  }
  counter=endLocInclusive; 
  //  cout<<"The buffer is "<<theValue;
}

void printTree(SectionClass *sec)
{
  cout<<"Section: "<<sec->Name<<endl;
  list<VarClass*>::iterator varIter=sec->VarList.begin();
  while (varIter!=sec->VarList.end()){
    cout<<"Variable: "<<(*varIter)->Name<<" ";
    VarASCIIClass *tempVar=(VarASCIIClass*)*varIter;
    printCharArray(tempVar->Value);
    //    cout<<tempVar->Value;
    varIter++;
  }
  list<SectionClass*>::iterator secIter=sec->SectionList.begin();
  while (secIter!=sec->SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    printTree(*secIter);
    secIter++;
  }
}


string currWord(Array <char,1> &buffer,
		int &counter,
		string &secName){
  string tempString="";
  int endLocExclusive=counter;
  bool inQuotes=false;
  while (isNumChar(buffer(endLocExclusive)) || inQuotes){
    tempString=tempString+buffer(endLocExclusive);
    if (buffer(endLocExclusive)=='\"' && buffer(endLocExclusive-1)!='\\'){
      inQuotes=!inQuotes;
    }
    endLocExclusive++;
  }
  ///We assume the section name doesn't contain bad characters inside quotes.
  if (tempString=="Section"){
    cerr<<"We have a section!!!"<<endl;
    while (!isNumChar(buffer(endLocExclusive))){
      endLocExclusive++;
    }
    secName="";
    while (isNumChar(buffer(endLocExclusive))){
      secName=secName+buffer(endLocExclusive);
      endLocExclusive++;
    }
    cout<<"The sec name is "<<secName<<endl;
  }
  counter=endLocExclusive-1;
  return tempString;
}

bool SpecialSectionClass<VarASCIIClass>::ReadFile(string fileName)
{


  SpecialSectionClass<VarASCIIClass> *sec;
  sec=this;
  SpecialSectionClass<VarASCIIClass> *oldSec=
    new SpecialSectionClass<VarASCIIClass>;
  //  oldSec=sec;
  //  VarASCIIClass *var = new VarASCIIClass;
  //  sec->Name = "Species";
  //  var->Name = "Mass";
  //  var->Value.resize(5);
  //  var->Value(0) = '3';
  //  var->Value(1) = '.';
  //  var->Value(2) = '1';
  //  var->Value(3) = '4';
  //  var->Value(4) = '2';
  //  sec->VarList.push_back(var);
  //  SectionList.push_back(sec);

  //  return true;
  string secName;
  int newCounter;
  Array<char,1> buffer;
  ReadWithoutComments(fileName,buffer);
  //  cout<<buffer;
  int braceLevel=0;
  bool inQuotes=false;
  sec->Name="all";

  for (int counter=0;counter<buffer.size();counter++){
    //    cout<<currWord(buffer,counter,secName)<<endl;
    if (buffer(counter)=='\"' && buffer(counter-1)!='\\'){
      inQuotes=!inQuotes;
    }
    else if (buffer(counter)=='S' &&
	     currWord(buffer,counter,secName)=="Section" &&
	     !inQuotes){
      SpecialSectionClass<VarASCIIClass> *newSec = 
	new SpecialSectionClass<VarASCIIClass>;
      newSec->Name=secName;
      newSec->parent=sec;
      sec->SectionList.push_back(newSec);
      sec=newSec;
    }
    else if (buffer(counter)=='}' &&
	     !inQuotes){
      //      sec=&((SpecialSectionClass<VarASCIIClass>)(*(sec->parent)))
      sec=(SpecialSectionClass<VarASCIIClass>*)(sec->parent);
    }
    else if(buffer(counter)=='='
	    && !inQuotes){
      string theName;
      Array <char,1> theValue;
      getVariable(buffer,counter,theName,theValue);
      VarASCIIClass *var = new VarASCIIClass;
      var->Name=theName;
      //      cout<<"Right here the value is "<<theValue;
      var->Value.resize(theValue.size());
      var->Value=theValue;
      sec->VarList.push_back(var);
    }
  }
  sec=oldSec;
  Iter =SectionList.begin();
}

//void PrintList(SectionClass *sec)
//{
  

//}






// int main()
// {
//   Array<char,1> myBuffer;
//   ReadWithoutComments("tryFile",myBuffer);
//   for (int counter=0;counter<myBuffer.size();counter++){
//     cout<<myBuffer(counter);
//   }

// }
 

void test()
{
  SpecialSectionClass<VarASCIIClass> section;
  SectionClass *newsec;
  section.ReadFile("myFile");
  cout<<"Here is "<<section.Name<<endl;
  cout<<"It is "<<section.SectionList.front()->Name<<"BBB"<<endl;
  bool found = section.FindSection("Species", newsec);
  if (found)
    {
      cerr << newsec->Name <<endl;      
      double x;
      //      found = newsec->ReadVar("Mass", x);
      if (found)
    	cerr << "x = " << x << endl;
    }
  
  

  cerr << found<<endl;
  printTree(&section);
  

}


main()
{
  test();
  exit(1);
}
