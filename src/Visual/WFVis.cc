#include "WFVis.h"
#include <GL/glut.h>
#include "ElementData.h"

WFVisualClass::WFVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  IsoAdjust  (0.01, 0.0, 1.0, 0.01, 0.1),
  BandAdjust (0.0, 0.0, 8.0, 1.0, 1.0),
  kAdjust (0.0, 0.0, 16.0, 1.0, 1.0),
  xPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  yPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  zPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  UpToDate(true),
  FileIsOpen(false),
  xPlane(WFIso),
  yPlane(WFIso),
  zPlane(WFIso),
  Export(*this),
  WFDisplay(MAG2),
  ResetIso(false)
{
  Glib::thread_init();
  WFIso.Dynamic = false;
  xPlane.Dynamic = false;
  yPlane.Dynamic = false;
  zPlane.Dynamic = false;
  IsoScale.set_adjustment(IsoAdjust);
  IsoScale.set_digits(2);
  IsoScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnIsoChange));
  IsoFrame.set_label("Density");
  IsoFrame.add(DensityBox);
  DensityBox.pack_start(rsLabel);
  DensityBox.pack_start(IsoScale);
  rsLabel.set_text("rs = ");
  IsoScale.property_width_request().set_value(75);

  BandScale.set_adjustment (BandAdjust);
  BandScale.set_digits(0);
  BandScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnBandChange));
  BandScale.property_width_request().set_value(75);
  BandFrame.set_label("Band");
  BandFrame.add(BandScale);

  kScale.set_adjustment (kAdjust);
  kScale.set_digits(0);
  kScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnkChange));
  kScale.property_width_request().set_value(75);
  kFrame.set_label("k point");
  kFrame.add(kScale);

  OrthoImage.set(FindFullPath("orthographic.png"));
  OrthoButton.set_icon_widget(OrthoImage);
  OrthoButton.set_label("Ortho");
  PerspectImage.set(FindFullPath("perspective.png"));
  PerspectButton.set_icon_widget(PerspectImage);
  PerspectButton.set_label("Perspect");

  ClipImage.set(FindFullPath("clipping.png"));
  ClipButton.set_icon_widget(ClipImage);
  ClipButton.set_label("Clip");

  IsoImage.set(FindFullPath("isoButton.png"));
  IsoButton.set_icon_widget(IsoImage);
  IsoButton.set_label("Isosurf");


  //////////////////////////////
  // Color plane widget setup //
  //////////////////////////////
  xPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  yPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  zPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  xPlaneScale.set_adjustment(xPlaneAdjust);
  yPlaneScale.set_adjustment(yPlaneAdjust);
  zPlaneScale.set_adjustment(zPlaneAdjust);
  xPlaneScale.property_width_request().set_value(75);
  yPlaneScale.property_width_request().set_value(75);
  zPlaneScale.property_width_request().set_value(75);
  xPlaneScale.set_digits(2);
  yPlaneScale.set_digits(2);
  zPlaneScale.set_digits(2);
  xPlaneButton.set_label("x Plane");
  yPlaneButton.set_label("y Plane");
  zPlaneButton.set_label("z Plane");
  xPlaneBox.pack_start (xPlaneButton, Gtk::PACK_SHRINK, 5);
  xPlaneBox.pack_start (xPlaneScale,  Gtk::PACK_SHRINK, 5);
  yPlaneBox.pack_start (yPlaneButton, Gtk::PACK_SHRINK, 5);
  yPlaneBox.pack_start (yPlaneScale,  Gtk::PACK_SHRINK, 5);
  zPlaneBox.pack_start (zPlaneButton, Gtk::PACK_SHRINK, 5);
  zPlaneBox.pack_start (zPlaneScale,  Gtk::PACK_SHRINK, 5);
  PlaneBox.pack_start(xPlaneBox, Gtk::PACK_SHRINK);
  PlaneBox.pack_start(yPlaneBox, Gtk::PACK_SHRINK);
  PlaneBox.pack_start(zPlaneBox, Gtk::PACK_SHRINK);
  PlaneFrame.add (PlaneBox);
  xPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 0));
  yPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 1));
  zPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 2));
  xPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 0));
  yPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 1));
  zPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 2));


  set_reallocate_redraws(true);
  PathVis.set_size_request(800, 800);
  ////////////////////
  // Setup tool bar //
  ////////////////////
  Gtk::RadioButtonGroup group = OrthoButton.get_group();
  PerspectButton.set_group(group);
  
  Tools.append(OrthoButton);
  Tools.append(PerspectButton);
  Tools.append(ClipButton);
  Tools.append(IsoButton);
  
  /////////////////
  // Setup menus //
  /////////////////
  Actions = Gtk::ActionGroup::create();
  Actions->add (Gtk::Action::create("MenuFile", "_File"));
  Actions->add (Gtk::Action::create("Open", "_Open"),
		sigc::mem_fun(*this, &WFVisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &WFVisualClass::OnExport));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &WFVisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &WFVisualClass::OnViewReset));
  CoordToggle = Gtk::ToggleAction::create("Axes", "Coordinate axes",
					  "Show coordinate axes", true);
  Actions->add (CoordToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnCoordToggle));
  SphereToggle = Gtk::ToggleAction::create("Nuclei", "Nuclei",
					  "Show nuclei", true);
  Actions->add (SphereToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnSphereToggle));

  Actions->add (Gtk::Action::create("MenuDisplay", "Display"));
  Mag2Radio = Gtk::RadioAction::create
    (DisplayGroup, "Mag2", "Magnitude squared",
     "Display square magnitude of WF");
  RealRadio = Gtk::RadioAction::create
    (DisplayGroup, "Real", "Real part", "Display real part of WF");
  ImagRadio = Gtk::RadioAction::create
    (DisplayGroup, "Imag", "Imag part", "Display imaginary part of WF");
  Actions->add
    (RealRadio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),REAL_PART));
  Actions->add
    (ImagRadio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),IMAG_PART));
  Actions->add
    (Mag2Radio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),MAG2));  

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
    "      <menuitem action='Nuclei'/>"
    "      <menuitem action='Axes'/>"
    "    </menu>"
    "    <menu action='MenuDisplay'>"
    "      <menuitem action='Real'/>"
    "      <menuitem action='Imag'/>"
    "      <menuitem action='Mag2'/>"
    "    </menu>"
    "  </menubar>"
    "  <toolbar  name='ToolBar'>"
    "    <toolitem action='Open'/>"
    "    <toolitem action='Quit'/>"
    "  </toolbar>"
    "</ui>";
  Manager = Gtk::UIManager::create();
  Manager->insert_action_group(Actions);
  add_accel_group (Manager->get_accel_group());
  Manager->add_ui_from_string (ui_info);

  /////////////////////
  // Connect signals //
  /////////////////////
  OrthoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnPerspectiveToggle));
  ClipButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnClipToggle));
  IsoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnIsoToggle));
  QuitButton.signal_clicked().connect
    (sigc::mem_fun(*this, &WFVisualClass::Quit));


  ////////////////////
  // Pack the boxes //
  ////////////////////
  ToolBox.pack_start(Tools);
  ToolBox.pack_start(PlaneFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(IsoFrame,  Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(BandFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(kFrame,    Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start(*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK,0);
  MainVBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(QuitButton, Gtk::PACK_SHRINK, 0);

  add (MainVBox);
  set_title ("wfvis++");
  show_all();

  PerspectButton.set_active(true);
  UpdateIso = true;
  UpdatePlane[0] = UpdatePlane[1] = UpdatePlane[2] = true;
}

void
WFVisualClass::OnViewReset()
{
  PathVis.View.Reset();
  PathVis.Invalidate();
}


void
WFVisualClass::OnOpen()
{

}

void
WFVisualClass::Quit()
{
  Gtk::Main::quit();
}

void
WFVisualClass::OnExport()
{
  Export.SetupWidgets();
  Export.show_all();
}





string 
WFVisualClass::FindFullPath(string filename)
{
  string fullName;

  fullName = filename;
  if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
    return fullName;
  else {
    fullName = PKG_DATA_DIR + filename;
    if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
      return fullName;
    else {
      cerr << "Cannot find \"" << filename << "\" anywhere.\n";
      return filename;
    }
  }
}



void
WFVisualClass::OnClipToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
}

void
WFVisualClass::OnIsoToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
}

void
WFVisualClass::OnPerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  DrawFrame();
  //  PathVis.Invalidate();
}


// void
// WFVisualClass::ClipSpheres(list<Vec3>& sphereList, double radius)
// {
//    list<Vec3>::iterator iter;
//    /// First, replicate spheres
//    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//      Vec3 &r = (*iter);
//      bool makeDisks = false;
//      if ((r[0]+radius) > 0.5*Box[0]) 
//        sphereList.push_front(Vec3(r[0]-Box[0], r[1], r[2]));
//       if ((r[0]-radius) < -0.5*Box[0]) 
// 	sphereList.push_front(Vec3(r[0]+Box[0], r[1], r[2]));
//     }
    
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       if ((r[1]+radius) > 0.5*Box[1])
// 	sphereList.push_front(Vec3(r[0], r[1]-Box[1], r[2]));
//       if ((r[1]-radius) < -0.5*Box[1])
// 	sphereList.push_front(Vec3(r[0], r[1]+Box[1], r[2]));
//     }
    
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       if ((r[2]+radius) > 0.5*Box[2])
// 	sphereList.push_front(Vec3(r[0], r[1], r[2]-Box[2]));
//       if ((r[2]-radius) < -0.5*Box[2])
// 	sphereList.push_front(Vec3(r[0], r[1], r[2]+Box[2]));
//     }
//     // Now make disks
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       for (int dim=0; dim<3; dim++) {
// 	if ((r[dim]+radius) > Box[dim]) {
// 	  double l = 0.5*Box[dim]-fabs(r[dim]);
// 	  double diskRad = sqrt(radius*radius-l*l);
// 	  DiskObject *disk1 = new DiskObject(offScreen);
// 	  DiskObject *disk2 = new DiskObject(offScreen);
// 	  disk1->SetRadius(diskRad);
// 	  disk2->SetRadius(diskRad);
// 	  disk1->SetAxis(dim);
// 	  disk2->SetAxis(-dim);
// 	  Vec3 r1, r2;
// 	  r1 = r; r2 = r;
// 	  r1[dim] =  0.5*Box[dim];
// 	  r2[dim] = -0.5*Box[dim];
// 	  disk1->SetPos(r1);
// 	  disk2->SetPos(r2);
// 	  PathVis.Objects.push_back(disk1);
// 	  PathVis.Objects.push_back(disk2);
//       }
// }

class AtomClass
{
public:
  Vec3 Pos;
  int Type;
  AtomClass(Vec3 pos, int type) {
    Pos = pos;
    Type = type;
  }
};

bool
WFVisualClass::DrawFrame(bool offScreen)
{
  bool clipping = ClipButton.get_active();
  for (int i=0; i<PathVis.Objects.size(); i++) {
    if (PathVis.Objects[i]->Dynamic) {
      if (dynamic_cast<SphereObject*> (PathVis.Objects[i]) != NULL)
	delete dynamic_cast<SphereObject*> (PathVis.Objects[i]);
      else if (dynamic_cast<DiskObject*> (PathVis.Objects[i]) != NULL)
	delete dynamic_cast<DiskObject*> (PathVis.Objects[i]);
      else if (dynamic_cast<Isosurface*> (PathVis.Objects[i]) != NULL)
	delete dynamic_cast<Isosurface*> (PathVis.Objects[i]);
      else if (dynamic_cast<BoxObject*> (PathVis.Objects[i]) != NULL)
	delete dynamic_cast<BoxObject*> (PathVis.Objects[i]);
    }

    //    delete (PathVis.Objects[i]);
  }
  PathVis.Objects.resize(0);

  if (CoordToggle->get_active()) {
    CoordObject *coord = new CoordObject;
    coord->Set (Box);
    PathVis.Objects.push_back(coord);
  }

  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  boxObject->Set (Box, clipping);
  PathVis.Objects.push_back(boxObject);
  
  
  if (SphereToggle->get_active()) {
    list<AtomClass>::iterator iter;
    list<AtomClass> sphereList;
    for (int ptcl=0; ptcl < AtomPos.extent(0); ptcl++) 
      sphereList.push_back(AtomClass(AtomPos(ptcl), AtomTypes(ptcl)));
    if (clipping) {
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = ElementData::GetRadius(type);
	bool makeDisks = false;
	if ((r[0]+radius) > 0.5*Box[0]) {
	  sphereList.push_front(AtomClass(Vec3(r[0]-Box[0],r[1], r[2]), type));
	  makeDisks = true;
	}
	if ((r[0]-radius) < -0.5*Box[0]) {
	  sphereList.push_front(AtomClass(Vec3(r[0]+Box[0],r[1], r[2]), type));
	  makeDisks = true;
	}
      }
      
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = ElementData::GetRadius(type);
	if ((r[1]+radius) > 0.5*Box[1])
	  sphereList.push_front(AtomClass(Vec3(r[0], r[1]-Box[1], r[2]),type));
	if ((r[1]-radius) < -0.5*Box[1])
	  sphereList.push_front(AtomClass(Vec3(r[0], r[1]+Box[1], r[2]),type));
      }
      
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = ElementData::GetRadius(type);
	if ((r[2]+radius) > 0.5*Box[2])
	  sphereList.push_front(AtomClass(Vec3(r[0], r[1], r[2]-Box[2]),type));
	if ((r[2]-radius) < -0.5*Box[2])
	  sphereList.push_front(AtomClass(Vec3(r[0], r[1], r[2]+Box[2]),type));
      }
      // Now make disks to close spheres
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = ElementData::GetRadius(type);
	
	for (int dim=0; dim<3; dim++) {
	  if ((r[dim]+radius) > 0.5*Box[dim]) {
	    double l = 0.5*Box[dim]-fabs(r[dim]);
	    double diskRad = sqrt(radius*radius-l*l);
	    DiskObject *disk1 = new DiskObject(offScreen);
	    DiskObject *disk2 = new DiskObject(offScreen);
	    Vec3 color = ElementData::GetColor (type);
	    disk1->SetRadius(diskRad);
	    disk2->SetRadius(diskRad);
	    disk1->SetAxis(2*dim);
	    disk2->SetAxis(2*dim+1);
	    disk1->SetColor(color);
	    disk2->SetColor(color);
	    Vec3 r1, r2;
	    r1 = r; r2 = r;
	    r1[dim] =  0.5*Box[dim];
	    r2[dim] = -0.5*Box[dim];
	    disk1->SetPos(r1);
	    disk2->SetPos(r2);
	    PathVis.Objects.push_back(disk1);
	    PathVis.Objects.push_back(disk2);
	  }
	}
      }
    }
    /// Add sphere objects from list
    for (iter=sphereList.begin(); iter!=sphereList.end(); iter++) {
      SphereObject *sphere = new SphereObject (offScreen);
      sphere->SetPos((*iter).Pos);
      Vec3 color = ElementData::GetColor (iter->Type);
      double radius = ElementData::GetRadius(iter->Type);
      sphere->SetColor(color);
      sphere->SetBox(Box);
      sphere->SetRadius(radius);
      PathVis.Objects.push_back(sphere);
    }
  }

  if (FileIsOpen) {
    int band, k;
    band = (int)round(BandScale.get_value());
    k    = (int)round(kScale.get_value());
    if ((CurrBand != band) || (Currk != k)) {
      CurrBand = band;
      Currk    = k;
      ReadWF (k, band);
      Xgrid.Init(-0.5*Box[0], 0.5*Box[0], WFData.extent(0));
      Ygrid.Init(-0.5*Box[1], 0.5*Box[1], WFData.extent(1));
      Zgrid.Init(-0.5*Box[2], 0.5*Box[2], WFData.extent(2));
      WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }
    if (ResetIso) {
      Xgrid.Init(-0.5*Box[0], 0.5*Box[0], WFData.extent(0));
      Ygrid.Init(-0.5*Box[1], 0.5*Box[1], WFData.extent(1));
      Zgrid.Init(-0.5*Box[2], 0.5*Box[2], WFData.extent(2));
      WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
      ResetIso = false;
    }
    if (UpdateIso) {
      if (WFDisplay == MAG2) {
	WFIso.SetColor (0.0, 0.8, 0.0);
	WFIso.SetIsoval(MaxVal*IsoAdjust.get_value());
      }
      else {
	vector<TinyVector<double,3> > colors;
	colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	WFIso.SetColor(colors);
	vector<double> vals;
	vals.push_back(+MaxVal*IsoAdjust.get_value());
	vals.push_back(-MaxVal*IsoAdjust.get_value());
	WFIso.SetIsoval(vals);
      }
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }

    if (UpdatePlane[0] && xPlaneButton.get_active())
      xPlane.SetPosition (0, xPlaneScale.get_value());
    if (xPlaneButton.get_active())
      PathVis.Objects.push_back(&xPlane);
    
    if (UpdatePlane[1] && yPlaneButton.get_active())
      yPlane.SetPosition (1, yPlaneScale.get_value());
    if (yPlaneButton.get_active())
      PathVis.Objects.push_back(&yPlane);
    
    if (UpdatePlane[2] && zPlaneButton.get_active())
      zPlane.SetPosition (2, zPlaneScale.get_value());
    if (zPlaneButton.get_active())
      PathVis.Objects.push_back(&zPlane);
    if (IsoButton.get_active()) 
      PathVis.Objects.push_back(&WFIso);
  }

  PathVis.Invalidate();
  UpToDate = true;
  UpdateIso = false;
  UpdatePlane[0] = false;
  UpdatePlane[1] = false;
  UpdatePlane[2] = false;
  return false;
}


bool 
IsDiag(Array<double,2> &lattice) 
{
  assert (lattice.extent(0) == 3);
  assert (lattice.extent(1) == 3);
  return ((fabs(lattice(0,1)) < 1.0e-16) &&
	  (fabs(lattice(1,0)) < 1.0e-16) &&
	  (fabs(lattice(0,2)) < 1.0e-16) &&
	  (fabs(lattice(2,0)) < 1.0e-16) &&
	  (fabs(lattice(1,2)) < 1.0e-16) &&
	  (fabs(lattice(2,1)) < 1.0e-16));
}


void
WFVisualClass::Read(string filename)
{
  if (FileIsOpen)
    Infile.CloseFile();
  assert (Infile.OpenFile(filename));
  FileIsOpen = true;
  
  /// Read lattice vectors
  assert (Infile.OpenSection("parameters"));
  Array<double,2> lattice;
  assert (Infile.ReadVar("lattice", lattice));
  if (!IsDiag(lattice)) { 
    cerr << "wfvis does not currently work with non-orthorhombic cells.\n";
    abort();
  }
  Box.Set (lattice(0,0), lattice(1,1), lattice(2,2));
  
  Infile.CloseSection(); // parameters

  /// Read ion positions
  assert (Infile.OpenSection("ions"));
  Array<double,2> pos;
  assert (Infile.ReadVar("pos", pos));
  AtomPos.resize(pos.extent(0));
  for (int i=0; i<pos.extent(0); i++)
    for (int j=0; j<3; j++)
      AtomPos(i)[j] = (pos(i,j) - 0.5*Box[j]);
  assert (Infile.ReadVar("atom_types", AtomTypes));
  Infile.CloseSection (); // "ions"

  /// Count k-points and bands
  assert (Infile.OpenSection("eigenstates"));
  Numk = Infile.CountSections ("twist");
  assert (Infile.OpenSection("twist", 0));
  NumBands = Infile.CountSections("band");
  Infile.CloseSection(); // "twist"
  Infile.CloseSection(); // "eigenstates"
  cerr << "Numk = " << Numk << "   NumBands = " << NumBands << endl;

  /// Read first wave function
  ReadWF(0,0);
  Xgrid.Init(-0.5*Box[0], 0.5*Box[0], WFData.extent(0));
  Ygrid.Init(-0.5*Box[1], 0.5*Box[1], WFData.extent(1));
  Zgrid.Init(-0.5*Box[2], 0.5*Box[2], WFData.extent(2));
  WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
  CurrBand = 0; 
  Currk = 0;
  BandAdjust.set_upper(NumBands-1.0);
  kAdjust.set_upper(Numk-1.0);

  double maxDim = max(max(lattice(0,0), lattice(1,1)), lattice(2,2));
  PathVis.View.SetDistance (1.2*maxDim);
  
  IsoButton.set_active(true);
  IsoButton.set_sensitive(true);
  IsoFrame.set_sensitive(true);
  DrawFrame();
}


	
void
WFVisualClass::OnIsoChange()
{
  double rho = IsoAdjust.get_value() * MaxVal;
  double rs = pow (3.0/(4.0*M_PI*rho),1.0/3.0);
  char rstext[100];
  snprintf (rstext, 100, "rs = %1.3f", rs);
  rsLabel.set_text(rstext);
  
  UpdateIso = true;
  DrawFrame();
}

void
WFVisualClass::OnBandChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  DrawFrame();
}

void
WFVisualClass::OnkChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  DrawFrame();
}

void
WFVisualClass::OnPlaneChange(int dir)
{
  if (dir==0 && xPlaneButton.get_active())
    UpdatePlane[0] = true;
  if (dir==1 && yPlaneButton.get_active())
    UpdatePlane[1] = true;
  if (dir==2 && zPlaneButton.get_active())
    UpdatePlane[2] = true;
  DrawFrame();
}


WFVisualClass::~WFVisualClass()
{

}

bool
WFVisualClass::ReadWF (int kpoint, int band)
{
  Array<double,4> wfdata;
  assert (Infile.OpenSection("eigenstates"));
  assert (Infile.OpenSection("twist", kpoint));
  assert (Infile.OpenSection("band", band));
  assert (Infile.ReadVar ("eigenvector", wfdata));
  Infile.CloseSection(); // "eigenstates"
  Infile.CloseSection(); // "twist"
  Infile.CloseSection(); // "band"
  WFData.resize(wfdata.extent(0), 
		wfdata.extent(1),
		wfdata.extent(2));
  int Nx = wfdata.extent(0);
  int Ny = wfdata.extent(1);
  int Nz = wfdata.extent(2);
  MaxVal = 0.0;
  for (int ix=0; ix<wfdata.extent(0); ix++)
    for (int iy=0; iy<wfdata.extent(1); iy++)
      for (int iz=0; iz<wfdata.extent(2); iz++) {
	// WF data is store from 0 to Lx, not -Lx/2 to Lx/2
	int jx = (ix/*+Nx/2*/)%(Nx-1);
	int jy = (iy/*+Ny/2*/)%(Ny-1);
	int jz = (iz/*+Nz/2*/)%(Nz-1);
	if (WFDisplay == MAG2) {
	  double rho = (wfdata(jx,jy,jz,0)*wfdata(jx,jy,jz,0) +
			wfdata(jx,jy,jz,1)*wfdata(jx,jy,jz,1));
	  WFData(ix,iy,iz) = rho;
	}
	else if (WFDisplay == REAL_PART)
	  WFData(ix,iy,iz) = wfdata(jx, jy, jz, 0);
	else if (WFDisplay == IMAG_PART)
	  WFData(ix,iy,iz) = wfdata(jx, jy, jz, 1);
	MaxVal = max(MaxVal, fabs(WFData(ix,iy,iz)));
      }
  return true;
}



void
WFVisualClass::OnCoordToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnSphereToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnDisplayRadio(WFDisplayType type)
{
  WFDisplayType newtype;
  if (RealRadio->get_active() && type==REAL_PART)
    newtype = REAL_PART;
  if (ImagRadio->get_active() && type==IMAG_PART)
    newtype = IMAG_PART;
  if (Mag2Radio->get_active() && type==MAG2)
    newtype = MAG2;

  if (WFDisplay != newtype) {
    WFDisplay = newtype;
    ReadWF (Currk, CurrBand);
    UpdatePlane[0] = UpdatePlane[1] = UpdatePlane[2] = true;
    UpdateIso = true;  ResetIso = true;
    DrawFrame ();
  }
}



int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);
  glutInit(&argc, argv);

  if (argc < 2) {
    cerr << "Usage:\n  wfvis++ myfile.h5\n";
    exit (1);
  }
  
  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  // Instantiate and run the application.
  WFVisualClass wfvisual;
  wfvisual.Read (argv[1]);
  kit.run(wfvisual);

  return 0;
}
