#include "GLObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>

using namespace std;

GLObject::GLObject() : Created(false), Visible(true)
{
 
};

GLObject::~GLObject()
{
  if (Created)
    glDeleteLists(ListNum, 1);
}

void GLObject::Start()
{
  if (Created)
    glDeleteLists(ListNum, 1);
  else
    ListNum = glGenLists(1);
  glNewList (ListNum, GL_COMPILE);
  Created = true;
}

void GLObject::End()
{
  glEndList();
}

void GLObject::Draw()
{
  if (Visible)
    glCallList (ListNum);
}

void GLObject::SetVisible (bool visible)
{
  Visible = visible;
}

bool GLObject::GetVisible ()
{
  return Visible;
}
