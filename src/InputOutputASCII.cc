#include "InputOutputASCII.h"

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
      theChar=='-')
    return true;
  return false;
}


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
  counter=endLocExclusive-1;
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


void SpecialInputSectionClass<VarASCIIClass>::printTree(InputSectionClass *sec )
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
  list<InputSectionClass*>::iterator secIter=sec->SectionList.begin();
  while (secIter!=sec->SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    printTree(*secIter);
    secIter++;
  }
}


void 
SpecialInputSectionClass<VarASCIIClass>::ReadWithoutComments(string fileName,
							     Array<char,1> 
							     &buffer)
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
}




bool SpecialInputSectionClass<VarASCIIClass>::ReadFile(string fileName)
{


  SpecialInputSectionClass<VarASCIIClass> *sec;
  sec=this;
  SpecialInputSectionClass<VarASCIIClass> *oldSec=
    new SpecialInputSectionClass<VarASCIIClass>;
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
      SpecialInputSectionClass<VarASCIIClass> *newSec = 
	new SpecialInputSectionClass<VarASCIIClass>;
      newSec->Name=secName;
      newSec->Iter=newSec->SectionList.begin();
      newSec->parent=sec;
      sec->SectionList.push_back(newSec);
      sec=newSec;
    }
    else if (buffer(counter)=='}' &&
	     !inQuotes){
      //      sec=&((SpecialSectionClass<VarASCIIClass>)(*(sec->parent)))
      sec=(SpecialInputSectionClass<VarASCIIClass>*)(sec->parent);
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



