#include "InputOutput.h"

void TestCommaSplit();
void TestStripWhiteSpace();

void TestXMLRead()
{
  IOSectionClass xmlFile;

  xmlFile.OpenFile ("Test1.xml");
  xmlFile.PrintTree();
  cerr << "Past PrintTree()\n";
}

main()
{
  TestCommaSplit();
  TestStripWhiteSpace();
  TestXMLRead();
}
