#ifndef WF_VIS_H
#define WF_VIS_H

#include "PathVis.h"
#include "BoxClass.h"
#include "Isosurface.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>


using namespace IO;

class WFVisualClass : public Gtk::Window
{
protected:
  //////////
  // Data //
  //////////
  Array<double,3> WFData;
  int CurrBand, Currk, NumBands, Numk;
  Array<Vec3,1> AtomPos;
  Array<int,1> AtomTypes;
  BoxClass Box;
  int Currentk, CurrentBand;

  /////////////
  // Widgets //
  /////////////
  Gtk::VBox MainVBox;
  Gtk::Button QuitButton;
  Gtk::HScale kScale, BandScale;
  Gtk::Adjustment kAdjust, BandAdjust;
  Gtk::Frame kFrame, BandFrame;
  
  Gtk::Toolbar Tools;
  Gtk::HBox ToolBox;
  Gtk::RadioToolButton OrthoButton, PerspectButton;
  Gtk::Image OrthoImage, PerspectImage;
  Gtk::ToggleToolButton ClipButton;
  Gtk::Image ClipImage;
  Gtk::ToggleToolButton IsoButton;
  Gtk::Image IsoImage;

  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;

  //////////////////////////////
  // Density isosurface stuff //
  //////////////////////////////
  IOSectionClass Infile;
  bool FileIsOpen;
  Array<double,3> RhoData;
  LinearGrid Xgrid, Ygrid, Zgrid;
  Gtk::VBox IsoBox, DensityBox;
  Gtk::Label rsLabel;
  Gtk::HScale IsoScale;
  Gtk::Adjustment IsoAdjust;
  Gtk::Frame IsoFrame;
  double FindMaxRho();
  double MaxRho;
  double FindMaxBandRho();
  double MaxBandRho;


  //////////////////////
  // Callback methods //
  //////////////////////
  void OnExport();
  void Quit();
  void OnSpeedChange();
  void OnIsoChange();
  void OnBandChange();
  void OnkChange();
  void OnPerspectiveToggle();
  void OnPlayToggle();
  void OnClipToggle();
  void OnIsoToggle();
  void OnViewReset();
  void OnOpen();
  bool UpToDate;

  ///////////////////
  // Other methods //
  ///////////////////
  string FindFullPath(string filename);
  bool   DrawFrame(bool offScreen=false);
  bool ReadWF (int kPoint, int band);
public:
  PathVisClass PathVis;

  ////////////////////
  // Public methods //
  ////////////////////
  void Read(string filename);
  
  WFVisualClass ();
  virtual ~WFVisualClass();
};

#endif
