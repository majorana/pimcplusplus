#ifndef VISUAL_H
#define VISUAL_H

#include "PathVis.h"
#include "../Common/IO/InputOutput.h"

class VisualClass : public Gtk::Window
{
protected:
  // signal handlers:
  void on_button_quit_clicked();

  // member widgets:
  Gtk::VBox m_VBox;
  Gtk::Button m_ButtonQuit;
  Array<double,4> Paths;

public:
  PathVisClass PathVis;

  void Read(string fileName);
  
  VisualClass();
  virtual ~VisualClass();

};

#endif
