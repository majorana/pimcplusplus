#ifndef EXPORT_H
#define EXPORT_H

#include <gtkmm.h>
#include <gtkglmm.h>

using namespace std;

class VisualClass;

class ExportClass : public Gtk::Window
{
private:
  VisualClass &Visual;
  Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap;
  Glib::RefPtr<Gdk::Pixmap>     GdkPixmap;
  Glib::RefPtr<Gdk::Drawable>   Draw;
  Glib::RefPtr<Gdk::GL::Config> GLConfig;
  Glib::RefPtr<Gdk::GL::Context> GLContext;
  int Width, Height;
  int AttribList;
  void InitGLStuff();
public:
  void Export (string filename);
  void ExportPOV (string filename);
  ExportClass (VisualClass &visual) : 
    Visual(visual), Width(2000), Height(2000)
  {
    
  };
};


#endif
