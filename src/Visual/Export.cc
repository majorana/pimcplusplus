#include "Export.h"
#include "Visual.h"
#include <GL/gl.h>
#include <GL/glu.h>


void ExportClass::InitGLStuff()
{
//   glShadeModel(GL_SMOOTH);
//   glEnable (GL_LIGHTING);
//   glEnable (GL_LINE_SMOOTH);
//   glEnable (GL_POLYGON_SMOOTH);

  static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
  static GLfloat light_specular[]= {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
  glEnable (GL_MULTISAMPLE);
  glEnable (GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  

  glViewport(0, 0, Width, Height);
  
}


void ExportClass::Export (string filename)
{
  GLConfig = Gdk::GL::Config::create (Gdk::GL::MODE_RGB   |
				      Gdk::GL::MODE_DEPTH |
				      Gdk::GL::MODE_SINGLE);
  if (!GLConfig) {
    cerr << "Cannot find a visual capable of OpenGL.\n";
    return;
  }

  GdkPixmap = Gdk::Pixmap::create (Visual.get_window(),
				   Width, Height, GLConfig->get_depth());
  
  Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap = 
    Gdk::GL::ext(GdkPixmap).set_gl_capability (GLConfig);

  GLContext = Gdk::GL::Context::create (GLPixmap, false);

  GLPixmap->make_current(GLContext);
  GLPixmap->gl_begin(GLContext);

  // Render here
  InitGLStuff();
  Visual.MakeFrame ((int)floor(Visual.FrameAdjust.get_value()));
  Visual.PathVis.GLRender();

  glFlush();
  GLPixmap->gl_end();
  GLPixmap->wait_gl();

  Glib::RefPtr<Gdk::Drawable> drawable =  GdkPixmap;
  Glib::RefPtr<Gdk::Drawable> drawable2 =  GdkPixmap;

  cerr << "drawable = " << drawable << endl;

  Glib::RefPtr<Gdk::Pixbuf> Pbuf = Gdk::Pixbuf::create
    (drawable, GdkPixmap->get_colormap(), 0, 0, 0, 0, Width, Height);

  if (drawable == drawable2)
    cerr << "no change\n";
						  
  Pbuf->save(filename, "png");

 //  GLPixmap = Gdk::GL::Pixmap::create 
//     (GLConfig, GdkPixmap, AttribList);
  
//  GLPixmap ->make_current(GLContext);


}



void
ExportClass::ExportPOV(string filename)
{
  Visual.PathVis.POVRender (filename);
}
