#ifndef MD_VIS_H
#define MD_VIS_H

#include "PathVis.h"
#include "BoxClass.h"
#include "Export.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>


using namespace IO;

class MDVisualClass : public Gtk::Window
{
protected:
  friend class ExportClass;
  //////////
  // Data //
  //////////
  Array<double,3> Trajectory;
  Array<double,1> Time;
  BoxClass Box;
  int CurrentFrame, PlayDirection, TimeoutDelay;

  /////////////
  // Widgets //
  /////////////
  Gtk::VBox MainVBox;
  Gtk::Button QuitButton;
  Gtk::HScale FrameScale;
  Gtk::Adjustment FrameAdjust;
  
  Gtk::Toolbar Tools;
  Gtk::HBox ToolBox;
  Gtk::RadioToolButton PlayButton, RevButton, PauseButton;
  Gtk::Image PlayImage, RevImage, PauseImage;
  Gtk::SeparatorToolItem ToolSep1;
  Gtk::RadioToolButton OrthoButton, PerspectButton;
  Gtk::Image OrthoImage, PerspectImage;

  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;

  //////////////////////
  // Callback methods //
  //////////////////////
  void OnExport();
  void Quit();
  void OnFrameChange();
  void OnPerspectiveToggle();
  void OnPlayToggle();
  void OnViewReset();
  void OnOpen();
  sigc::connection TimeoutConnection;
  /// The timeout callback is used for animation
  bool OnTimeout();

  ///////////////////
  // Other methods //
  ///////////////////
  string FindFullPath(string filename);
  void   DrawFrame();
public:
  PathVisClass PathVis;

  ////////////////////
  // Public methods //
  ////////////////////
  void Read(string filename);
  
  MDVisualClass ();
  virtual ~MDVisualClass();
};

#endif
