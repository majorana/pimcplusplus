#include <Python.h>
#include "InputOutput.h"



extern "C" PyObject* IOSection_New(PyObject *self, PyObject *args)
{
  IOSectionClass *newSect=new IOSectionClass();
  return Py_BuildValue("i", newSect);
}

extern "C" PyObject* IOSection_ReadVarInt(PyObject *self, PyObject *args)
{
  char* varName;
  void* IOSectionPtr;
  if (!PyArg_ParseTuple (args, "is", &IOSectionPtr,&varName))
      return NULL;
  else {
    int val;
    ((IOSectionClass*)IOSectionPtr)->ReadVar(varName,val);
    return Py_BuildValue("i",val);
  }    
}


extern "C" PyObject*
IOSection_OpenFile (PyObject *self, PyObject *args)
{
  char *fileName;
  void *IOSectionPtr;

  if (!PyArg_ParseTuple (args, "is",&IOSectionPtr,&fileName))
    return NULL;
  else {
    bool success = ((IOSectionClass*)IOSectionPtr)->OpenFile(fileName);
    return Py_BuildValue("i",(int)success);
  }

}
