#include "InputOutput.h"

#include <libxml/parser.h>




void startElementWrapper(void *parserPtr, const xmlChar *name,
			 const xmlChar **attr)
{
  list<XMLattribute> attributes;
  string nameStr = (const char *)name;
  int i=0;
  XMLattribute myAttr;
  while (attr[2*i] != NULL) {
    myAttr.Name = (const char *)(attr[i]);
    myAttr.Value = (const char *)(attr[i+1]);
    attributes.push_back(myAttr);
  }
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.startElement(nameStr, attributes);
}


void endElementWrapper(void *parserPtr, const xmlChar *name)
{
  string nameStr = (const char *)name;
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.endElement(nameStr);
}


void charactersWrapper(void *parserPtr, const xmlChar *ch, int len)
{
  string str((const char *)ch, len);
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.characters(str);
}



void XMLparserClass::SetHandler()
{
  handler.internalSubset        = NULL;
  handler.isStandalone          = NULL;
  handler.hasInternalSubset     = NULL;
  handler.hasExternalSubset     = NULL;
  handler.resolveEntity         = NULL;
  handler.getEntity             = NULL;
  handler.entityDecl            = NULL;
  handler.notationDecl          = NULL;
  handler.attributeDecl         = NULL;
  handler.elementDecl           = NULL;
  handler.unparsedEntityDecl    = NULL;
  handler.setDocumentLocator    = NULL;
  handler.startDocument         = NULL;
  handler.endDocument           = NULL;
  handler.startElement          = startElementWrapper;
  handler.endElement            = endElementWrapper;
  handler.reference             = NULL;
  handler.characters            = charactersWrapper;
  handler.ignorableWhitespace   = NULL;
  handler.processingInstruction = NULL;
  handler.comment               = NULL;
  handler.warning               = NULL;
  handler.error                 = NULL;
  handler.fatalError            = NULL;
}


bool XMLparserClass::ParseFile(string fileName, IOTreeXMLClass *rootNode)
{
  // set up the handler functions
  SetHandler();
  // store our filename
  currTree->FileName = fileName;
  
  int retval = xmlSAXUserParseFile(&handler, this, fileName.c_str());
  // don't know if this is right.
  return (retval == 0);
}


void XMLparserClass::startElement(string &name, list<XMLattribute> &attributes)
{
  cerr << "Element name:  " << name << endl;
  cerr << "Attributes:  " << endl;
  list<XMLattribute>::iterator iter = attributes.begin();
  while (iter != attributes.end()) {
    cerr << iter->Name << "=\"" << iter->Value << "\"\n";
    iter++;
  }
}

void XMLparserClass::endElement(string &name)
{
  cerr << "End element " << name << endl;
}

void XMLparserClass::characters(string &newChars)
{
  cerr << newChars;
}



bool IOTreeXMLClass::OpenFile(string fileName, string myName,
			      IOTreeClass *parent)
{
  Parent = parent;
  Name = myName;
  XMLparserClass parser;
  return parser.ParseFile (fileName, this);
}




/// Simply prints 3*num spaces
inline void XMLPrintIndent(int num)
{
  for (int counter=0;counter<num*2;counter++)
    cout<<' ';
}

/// Simply prints 3*num spaces
inline void XMLPrintIndent(int num,ofstream &outFile)
{
  for (int counter=0;counter<num*2;counter++)
    outFile<<' ';
}


/// Prints an indented hierarchy of sections and variable names to
/// cout. 
void IOTreeXMLClass::PrintTree(int indentNum)
{
  XMLPrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    XMLPrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

/// Calls PrintTree(0)
void IOTreeXMLClass::PrintTree()
{
  PrintTree(0);
}





VarXMLClass *NewXMLVar (AtomicType newType, int ndim,
			    Array<int,1> dims)
{
  if (ndim == 0) {
    if (newType == DOUBLE_TYPE)
      return new VarXMLdouble0Class;
    else if (newType == INT_TYPE)
      return new VarXMLint0Class;
    else if (newType == STRING_TYPE)
      return new VarXMLstring0Class;
    else if (newType == BOOL_TYPE)
      return new VarXMLbool0Class;
  }
  else if (ndim == 1) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble1Class *newVar = new VarXMLdouble1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint1Class *newVar = new VarXMLint1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring1Class *newVar = new VarXMLstring1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool1Class *newVar = new VarXMLbool1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
  }
  else if (ndim == 2) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble2Class *newVar = new VarXMLdouble2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint2Class *newVar = new VarXMLint2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring2Class *newVar = new VarXMLstring2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool2Class *newVar = new VarXMLbool2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
  }  
  else if (ndim == 3) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble3Class *newVar = new VarXMLdouble3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint3Class *newVar = new VarXMLint3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring3Class *newVar = new VarXMLstring3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool3Class *newVar = new VarXMLbool3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
  }  
}








IOTreeClass* IOTreeXMLClass::NewSection(string name)
{
  IOTreeClass* tempSection=new IOTreeXMLClass();
  tempSection->Name=name;
  tempSection->Parent=this;
  tempSection->MyNumber=CurrSecNum;
  CurrSecNum++;
  SectionList.push_back(tempSection);
  MarkModified();
  return tempSection;
}

void IOTreeXMLClass::IncludeSection(IOTreeClass *newSection)
{
  newSection->MyNumber=CurrSecNum++;
  SectionList.push_back(newSection);
  MarkModified();
}



bool IOTreeXMLClass::NewFile (string fileName,
			      string mySectionName,
			      IOTreeClass *parent)
{
  FileName=fileName;
  Parent=parent;
  Name=mySectionName;
}


void IOTreeXMLClass::WriteSection(ofstream &outFile,int indentNum)
{
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    XMLPrintIndent(indentNum,outFile);
    ((VarXMLClass*)(*varIter))->Print(outFile);    
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    if ((*secIter)->FileName==""){
      XMLPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<")\n";
      XMLPrintIndent(indentNum,outFile);
      outFile<<"{\n";
      ((IOTreeXMLClass*)(*secIter))->WriteSection(outFile,indentNum+1);
      XMLPrintIndent(indentNum,outFile);
      outFile<<"}\n\n";
    }
    else {
      XMLPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<", \"";
      outFile<<(*secIter)->FileName<<"\");"<<endl;
      (*secIter)->FlushFile();
    }
    secIter++;
  }
}



void IOTreeXMLClass::FlushFile()
{
  ofstream outfile;
  if ((FileName!="") && IsModified){
    outfile.open(FileName.c_str());
    WriteSection(outfile,0);
  }

  list<IOTreeClass*>::iterator iter = SectionList.begin();
  while (iter != SectionList.end()) {
    (*iter)->FlushFile();
    iter++;
  }
}

/// CloseFile recursively destroys the tree of data we constructed.
void IOTreeXMLClass::CloseFile()
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
}



//////////////////////////////////////////////////////////////////////
//                             ReadInto's                           //
//////////////////////////////////////////////////////////////////////

bool VarXMLClass::ReadInto (double &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (int &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (string &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (bool &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<double,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<double,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<double,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<int,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<int,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<int,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<string,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<string,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<string,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<bool,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<bool,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (Array<bool,3> &val)
{ ComplainReadInto(); return false; }

bool VarXMLClass::Append (double val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (int val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (string val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (bool val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<double,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<double,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<int,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<int,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<string,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<string,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<bool,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (Array<bool,2> &val)
{ ComplainAppend(); return false; }


bool VarXMLdouble0Class::ReadInto (double &val)
{ val = Value; return true;}
bool VarXMLint0Class::ReadInto (int &val)
{  val = Value; return true; }
bool VarXMLstring0Class::ReadInto (string &val)
{  val = Value; return true; }
bool VarXMLbool0Class::ReadInto (bool &val)
{ val = Value; return true; }
bool VarXMLdouble1Class::ReadInto (Array<double,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLdouble2Class::ReadInto (Array<double,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLdouble3Class::ReadInto (Array<double,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); val = Value; return true; }
bool VarXMLint1Class::ReadInto (Array<int,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLint2Class::ReadInto (Array<int,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLint3Class::ReadInto (Array<int,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarXMLstring1Class::ReadInto (Array<string,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLstring2Class::ReadInto (Array<string,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLstring3Class::ReadInto (Array<string,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarXMLbool1Class::ReadInto (Array<bool,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLbool2Class::ReadInto (Array<bool,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLbool3Class::ReadInto (Array<bool,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }


/*****************************************************************
 *                         Variable  Appends                     *
 *****************************************************************/
bool VarXMLdouble1Class::Append (double val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLint1Class::Append (int val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLstring1Class::Append (string val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLbool1Class::Append (bool val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLdouble2Class::Append (Array<double,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,Range::all()) = val;
  return(true);
}
bool VarXMLdouble3Class::Append (Array<double,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,Range::all(),Range::all()) = val;
  return(true);
}
bool VarXMLint2Class::Append (Array<int,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,Range::all()) = val;
  return(true);
}
bool VarXMLint3Class::Append (Array<int,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,Range::all(),Range::all()) = val;
  return(true);
}
bool VarXMLstring2Class::Append (Array<string,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,Range::all()) = val;
  return(true);
}
bool VarXMLstring3Class::Append (Array<string,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,Range::all(),Range::all()) = val;
  return(true);
}
bool VarXMLbool2Class::Append (Array<bool,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,Range::all()) = val;
  return(true);
}
bool VarXMLbool3Class::Append (Array<bool,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,Range::all(),Range::all()) = val;
  return(true);
}

void VarXMLdouble0Class::Print (ofstream &outFile)
{ outFile << "double " << Name << " = " << Value << ";\n";}
void VarXMLint0Class::Print (ofstream &outFile)
{ outFile << "int " << Name << " = " << Value << ";\n";} 
void VarXMLstring0Class::Print (ofstream &outFile)
{ outFile << "string " << Name << " = " << "\"" << Value << "\";\n";} 
void VarXMLbool0Class::Print (ofstream &outFile)
{ outFile << "bool " << Name <<" = "<< (Value ? "true;\n" : "false;\n"); }
void VarXMLdouble1Class::Print (ofstream &outFile)
{ 
  outFile << "Array<double,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarXMLdouble2Class::Print (ofstream &outFile)
{
  outFile << "Array<double,2> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLdouble3Class::Print (ofstream &outFile)
{
  outFile << "Array<double,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarXMLint1Class::Print (ofstream &outFile)
{
  outFile << "Array<int,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarXMLint2Class::Print (ofstream &outFile)
{
  outFile << "Array<int,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLint3Class::Print (ofstream &outFile)
{
  outFile << "Array<int,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarXMLstring1Class::Print (ofstream &outFile)
{
  outFile << "Array<string,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << "\"" << Value(i) << "\"" << ", ";
  outFile << "\"" << Value(Value.extent(0)-1) << "\"];\n";
}
void VarXMLstring2Class::Print (ofstream &outFile)
{
  outFile << "Array<string,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << "\"" << Value(i,j) << "\", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLstring3Class::Print (ofstream &outFile)
{
  outFile << "Array<string,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << "\"" << Value(i,j,k) << "\", ";
  outFile << "\"" 
	  << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "\"];\n";
}
void VarXMLbool1Class::Print (ofstream &outFile)
{
  outFile << "Array<bool,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << (Value(i) ? "true" : "false") << ", ";
  outFile << (Value(Value.extent(0)-1) ? "true" : "false") << "];\n";
}
void VarXMLbool2Class::Print (ofstream &outFile)
{
  outFile << "Array<bool,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << (Value(i,j) ? "true" : "false") << ", ";
  outFile << 
    (Value(Value.extent(0)-1,Value.extent(1)-1) ? "true" : "false") << "];\n";
}
void VarXMLbool3Class::Print (ofstream &outFile)
{
  outFile << "Array<bool,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << (Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) ?
	      "true" : "false") << "];\n";
}




//////////////////////////////////////////////////////////////////////
//                             Append's                             //
//////////////////////////////////////////////////////////////////////


void IOTreeXMLClass::WriteVar(string name, double val)
{
  VarXMLdouble0Class *newVar = new VarXMLdouble0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<double,1> &val)
{
  VarXMLdouble1Class *newVar = new VarXMLdouble1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<double,2> &val)
{
  VarXMLdouble2Class *newVar = new VarXMLdouble2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<double,3> &val)
{
  VarXMLdouble3Class *newVar = new VarXMLdouble3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeXMLClass::WriteVar(string name, int val)
{
  VarXMLint0Class *newVar = new VarXMLint0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<int,1> &val)
{
  VarXMLint1Class *newVar = new VarXMLint1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<int,2> &val)
{
  VarXMLint2Class *newVar = new VarXMLint2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<int,3> &val)
{
  VarXMLint3Class *newVar = new VarXMLint3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeXMLClass::WriteVar(string name, string val)
{
  VarXMLstring0Class *newVar = new VarXMLstring0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<string,1> &val)
{
  VarXMLstring1Class *newVar = new VarXMLstring1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<string,2> &val)
{
  VarXMLstring2Class *newVar = new VarXMLstring2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<string,3> &val)
{
  VarXMLstring3Class *newVar = new VarXMLstring3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
 
void IOTreeXMLClass::WriteVar(string name, bool val)
{
  VarXMLbool0Class *newVar = new VarXMLbool0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<bool,1> &val)
{
  VarXMLbool1Class *newVar = new VarXMLbool1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<bool,2> &val)
{
  VarXMLbool2Class *newVar = new VarXMLbool2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, Array<bool,3> &val)
{
  VarXMLbool3Class *newVar = new VarXMLbool3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
