#include "InputOutputASCII.h"
#include "InputOutputHDF5.h"








bool checkPair(Array<char,1> &buffer,int counter,char* toSee)
{
  if (counter+1>=buffer.size()){
    return false;
  }
  if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
    return true;
  }
  else return false;

}
bool isNumCharP(char theChar)
{
  if ((theChar>='A' && theChar<='Z') ||
      (theChar>='a' && theChar<='z') ||
      (theChar>='0' && theChar<='9') ||
      theChar=='_' ||
      theChar=='\"'||
      theChar=='-' ||
      theChar=='.')
    return true;
  return false;
}

//Gets the current word and puts the counter at the end of it if it is not
//section. Otherwise puts the counter at the end after the section name.
string currWord(Array<char,1> &buffer,int &counter,string &secName)

{
  string tempString="";
  int endLocExclusive=counter;
  bool inQuotes=false;
  while (isNumCharP(buffer(endLocExclusive)) || inQuotes){
    tempString=tempString+buffer(endLocExclusive);
    if (buffer(endLocExclusive)=='\"' && buffer(endLocExclusive-1)!='\\'){
      inQuotes=!inQuotes;
    }
    endLocExclusive++;
  }
  ///We assume the section name doesn't contain bad characters inside quotes.
  if (tempString=="Section"){
    cerr<<"We have a section!!!"<<endl;
    while (!isNumCharP(buffer(endLocExclusive))){
      endLocExclusive++;
    }
    secName="";
    while (isNumCharP(buffer(endLocExclusive))){
      secName=secName+buffer(endLocExclusive);
      endLocExclusive++;
    }
    cout<<"The sec name is "<<secName<<endl;
  }
  counter=endLocExclusive;
  return tempString;
}


void printCharArray(Array <char,1> theArray)
{
  for (int counter=0;counter<theArray.size();counter++){
    cout<<theArray(counter);
  }
  cout<<endl;
}

void getVariable(Array <char,1> &buffer,int &counter,string &theName,
		 Array <char,1> &theValue)
{
  theName="";
  int eqLocation=counter;
  while ((!(isNumCharP(buffer(counter)))) &&
	 counter>0){
    counter--;
  }

  int endLocInclusive=counter;
  while ((isNumCharP(buffer(counter)))
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
  while (!isNumCharP(buffer(beginLocInclusive))){
    beginLocInclusive++;
  }
  endLocInclusive=beginLocInclusive;
  while (isNumCharP(buffer(endLocInclusive)) ||  inQuotes){
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

void InputSectionASCIIClass::PrintTree(int num)
{
  cerr<<"This shouldn't be called yet!"<<endl;
}

void InputSectionASCIIClass::PrintTree()
{

  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    cout<<"Variable: "<<(*varIter)->Name<<" ";
    VarASCIIClass *tempVar=(VarASCIIClass*)*varIter;
    printCharArray(tempVar->Value);
    //    cout<<tempVar->Value;
    varIter++;
  }
  list<InputSectionClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree();
    secIter++;
  }
}

void InputSectionASCIIClass::CloseFile()
{
  return;

}
void 
InputSectionASCIIClass::ReadWithoutComments(string fileName)
							    
{
  ifstream infile;
  infile.open(fileName.c_str());
  Array<char,1> tmpBuffer;
  int counter=0;
  pos=0;
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
      while ( counter<tmpBuffer.size() && tmpBuffer(counter) != '\"'){
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
  infile2.close();
}



bool isWhiteSpace(char testChar)
{
  if (testChar==' ' || 
      testChar=='\n' ||
      testChar=='\t' ||
      testChar=='\r')
    return true;
  return false;
}

//Leaves at the first character that's not white space
void InputSectionASCIIClass::SkipWhiteSpace()
{
  while (isWhiteSpace(buffer(pos))){
    pos++;
  }
  return;
}

bool isNameChar(char testChar)
{
  return false ||
    (testChar>='0' && testChar<='9') ||
    (testChar>='a' && testChar<='z') ||
    (testChar>='A' && testChar<='Z') ||
    (testChar=='_');   
}


///Leaves pos at the first point outside the word
string InputSectionASCIIClass::ReadCurrWord()
{
  SkipWhiteSpace();
  string returnString="";
  while (pos<buffer.size() &&
	 isNameChar(buffer(pos))){
    returnString=returnString+buffer(pos);
    pos++;
  }
  return returnString;
}



string InputSectionASCIIClass::ReadValueString()
{
  bool inQuotes=false;
  string returnString;
  while (pos<buffer.size() &&
	 (inQuotes ||
	 buffer(pos)!=';')){
    if (buffer(pos)=='\"'){
      inQuotes=!inQuotes;
    }
    returnString=returnString+buffer(pos);
    pos++;
  }
  assert(pos>=buffer.size());
  return returnString;
}


AtomicType NameToType (string typeName)
{
  switch (typeName) {
  case "int":
    return (INT_TYPE);
  case "double":
    return (DOUBLE_TYPE);
  case "string":
    return (STRING_TYPE);
  case: "bool":
    return (BOOL_TYPE);
  default:
    return (NOT_ATOMIC);
  }
}


class ASCIIVarTypeClass
{
  virtual StringToValue


}

void ReadVal(string buf, int &A)
{

}

template <class T>
void ReadVal (string buf, Array<T,1> &destArray)
{


}

ReadVal (buf, Array<Array<int,1>,1> myarray)
{
  
  

}

InputSectionASCIIClass::ParseArray()
{
  SkipWhiteSpace();
  assert(buffer(pos)=='<');
  string typeName=GetCurrWord();
  AtomicType myType=NameToType(typeName);
  assert (myType!=NOT_ATOMIC);
  SkipWhiteSpace();
  assert(buffer(pos)==',');
  SkipWhiteSpace();
  string intString=GetCurrWord();
  char *ptr;
  int dim=strtol(intString.c_str(),&ptr,10);
  assert  (*ptr==NULL);
  SkipWhiteSpace();
  assert(buffer(pos)=='>');
  string myName=GetCurrWord();
  SkipWhiteSpace();
  assert(buffer(pos)=='=');
  SkipWhiteSpace();
  string buf=ReadValueString();
  
  ASCIIVarClass *newVar=new ASCIIVarClass();
  newVar->Dim=dim;
  newVar->Type=myType;
  switch(myType){
  case INT_TYPE:
    if (dim==1){
      Array<int,1> *tempArray=new Array<int,1>();
      
      ReadVal(
    }
    else if (dim==2){
    }
    else if (dim==3){
    }

  
}


InputSectionASCIIClass::ReadVariableInfo(string theType)
{
  if (NameToType(theType)==NOT_ATOMIC){
    assert(theType=="Array");
    ParseArray()
      SkipWhiteSpace();
    assert(buffer(pos)=='<');
    SkipWhiteSpace();
    
    
      
  }
  else {
  }
    


}


bool InputSectionASCIIClass::ReadSection(string sectionName,
					 InputSectionClass* parent,
					 bool findBrace)
{
  Name=sectionName;
  Parent=parent;
  while (!done){
    string myCurrWord=ReadCurrWord();
    SkipWhiteSpace();
    if (myCurrWord=="Section" &&
	buffer(pos)=='('){
      pos++;
      string theName=ReadCurrWord();
      SkipWhiteSpace();
      assert(buffer(pos)==')');
      pos++;
      SkipWhiteSpace();
      assert(buffer(pos)=='{');
      pos++;
      InputSectionClass *newSec= new InputSectionASCIIClass();
      SectionList.push_back(newSec);
      (*newSec)->ReadSection(theName,this);
      
    }
    else {
      ReadVariableInfo(myCurrWord);
    }
}
    
    





}

bool InputSectionASCIIClass::OpenFile(string fileName,
				      string sectionName,
				      InputSectionClass* parent)
{
  
  ReadWithoutComments(fileName);
  ReadSection(sectionName,parent,false);
  return true;


}


////Section Name PArt not implemented yet
bool InputSectionASCIIClass::OpenFile(string fileName, 
				      string sectionName, 
				      InputSectionClass *parent)
{
  InputSectionASCIIClass *sec;
  sec=this;
  InputSectionASCIIClass *oldSec=
    new InputSectionASCIIClass;
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
  sec->Parent=parent;
  for (int counter=0;counter<buffer.size();counter++){
    //    cout<<currWord(buffer,counter,secName)<<endl;
    if (buffer(counter)=='\"' && buffer(counter-1)!='\\'){
      inQuotes=!inQuotes;
    }
    else if (buffer(counter)=='S' &&
	     currWord(buffer,counter,secName)=="Section" &&
	     !inQuotes){
      InputSectionASCIIClass *newSec = 
	new InputSectionASCIIClass;
      newSec->Name=secName;
      newSec->Iter=newSec->SectionList.begin();
      newSec->Parent=sec;
      sec->SectionList.push_back(newSec);
      sec=newSec;
    }
    else if (buffer(counter)=='}' &&
	     !inQuotes){
      //      sec=&((SpecialSectionClass<VarASCIIClass>)(*(sec->Parent)))
      sec=(InputSectionASCIIClass*)(sec->Parent);
    }
    else if (buffer(counter)=='d' &&
	     currWord(buffer,counter,secName)=="double" &&
	     !inQuotes){
      string dummyString;
      VarASCIIClass *var =new VarASCIIClass();
      while (!isNumCharP(buffer(counter))){
	counter++;
      }
      var->Name=currWord(buffer,counter,dummyString);
      cerr<<"My name is "<<var->Name<<endl;
      while (!isNumCharP(buffer(counter))){
	counter++;
      }

      string tempValue;
      tempValue=currWord(buffer,counter,dummyString);
      cerr<<"The string of my double is "<<tempValue<<endl;
      var->PtrValue=new double(atof(tempValue.c_str()));
      sec->VarList.push_back(var);
    } 
    else if (buffer(counter)=='i' &&
	     currWord(buffer,counter,secName)=="int" &&
	     !inQuotes){
      string dummyString;
      VarASCIIClass *var =new VarASCIIClass();
      while (!isNumCharP(buffer(counter))){
	counter++;
      }
      var->Name=currWord(buffer,counter,dummyString);
      cerr<<"My name is "<<var->Name<<endl;
      while (!isNumCharP(buffer(counter))){
	counter++;
      }

      string tempValue;
      tempValue=currWord(buffer,counter,dummyString);
      cerr<<"The string of my int is "<<tempValue<<endl;
      var->PtrValue=new double(atoi(tempValue.c_str()));
      sec->VarList.push_back(var);
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
  //  sec=oldSec;
  Iter =SectionList.begin();
  return true;
}



