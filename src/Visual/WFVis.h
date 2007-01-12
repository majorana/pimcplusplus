#ifndef WF_VIS_H
#define WF_VIS_H

#include "PathVis.h"
#include "BoxClass.h"
#include "CoordObject.h"
#include "Isosurface.h"
#include "PlaneObject.h"
#include "WFExport.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>


using namespace IO;

typedef enum {REAL_PART, IMAG_PART, MAG2} WFDisplayType;

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
  WFDisplayType WFDisplay;

  /////////////
  // Widgets //
  /////////////
  Gtk::VBox MainVBox;
  Gtk::Button QuitButton;
  Gtk::HScale kScale, BandScale;
  Gtk::Adjustment kAdjust, BandAdjust;
  Gtk::Frame kFrame, BandFrame;
  Glib::RefPtr<Gtk::ToggleAction> CoordToggle, SphereToggle;
  Gtk::RadioButtonGroup DisplayGroup;
  Glib::RefPtr<Gtk::RadioAction> RealRadio, ImagRadio, Mag2Radio;
  
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
  Isosurface WFIso;
  LinearGrid Xgrid, Ygrid, Zgrid;
  Gtk::VBox IsoBox, DensityBox;
  Gtk::Label rsLabel;
  Gtk::HScale IsoScale;
  Gtk::Adjustment IsoAdjust;
  Gtk::Frame IsoFrame;
  double FindMaxVal();
  double MaxVal;
  double FindMaxBandVal();
  double MaxBandVal;

  ///////////////////////
  // Color plane stuff //
  ///////////////////////
  PlaneObject xPlane, yPlane, zPlane;
  Gtk::Frame PlaneFrame;
  Gtk::VBox PlaneBox;
  Gtk::Adjustment xPlaneAdjust, yPlaneAdjust, zPlaneAdjust;
  Gtk::HScale xPlaneScale, yPlaneScale, zPlaneScale;
  Gtk::CheckButton xPlaneButton, yPlaneButton, zPlaneButton;
  Gtk::HBox xPlaneBox, yPlaneBox, zPlaneBox;

  /////////////////
  // State flags //
  /////////////////
  bool UpdateIso, ResetIso;
  bool UpdatePlane[3];
  
  /////////////////////////////////////
  // Saving and opening viewer state //
  /////////////////////////////////////
  bool WriteState(string fname);
  void OnOpenState();
  void OnSaveState();
  Gtk::FileChooserDialog OpenStateChooser, SaveStateChooser;


  ////////////
  // Export //
  ////////////
  WFExportClass Export;

  //////////////////////
  // Callback methods //
  //////////////////////
  void OnExport();
  void Quit();
  void OnSpeedChange();
  void OnIsoChange();
  void OnBandChange();
  void OnkChange();
  void OnPlaneChange(int dir);
  void OnPerspectiveToggle();
  void OnPlayToggle();
  void OnClipToggle();
  void OnIsoToggle();
  void OnViewReset();
  void OnCoordToggle();
  void OnSphereToggle();
  void OnOpen();
  void OnDisplayRadio(WFDisplayType type);
  bool UpToDate;

  ///////////////////
  // Other methods //
  ///////////////////
  string FindFullPath(string filename);
  bool ReadWF (int kPoint, int band);
public:
  PathVisClass PathVis;

  ////////////////////
  // Public methods //
  ////////////////////
  void Read(string filename);
  bool ReadState (string fname);
  bool DrawFrame(bool offScreen=false);

  
  WFVisualClass ();
  virtual ~WFVisualClass();
};

#endif
