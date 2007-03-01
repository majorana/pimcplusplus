#include "WFVis.h"
#include <GL/glut.h>
#include "ElementData.h"
#include "ParseCommand.h"



WFVisualClass::WFVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  IsoAdjust  (0.01, 0.0, 1.0, 0.01, 0.1),
  BandAdjust (0.0, 0.0, 8.0, 1.0, 1.0),
  kAdjust (0.0, 0.0, 16.0, 1.0, 1.0),
  xPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  yPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  zPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  RadiusAdjust(0.4, 0.0, 1.0, 0.01, 0.1),
  UpToDate(true),
  FileIsOpen(false),
  xPlane(WFIso),
  yPlane(WFIso),
  zPlane(WFIso),
  Export(*this),
  WFDisplay(MAG2),
  ResetIso(false),
  SaveStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
  OpenStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_OPEN),
  UpdateIsoType(false), 
  UpdateIsoVal(false),
  DoShift (false),
  Shift (0.0, 0.0, 0.0),
  CMap(BLUE_WHITE_RED)
{
  //Glib::thread_init();
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

  RadiusScale.set_adjustment(RadiusAdjust);
  RadiusScale.set_digits(2);
  RadiusScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnRadiusChange));
  RadiusScale.property_width_request().set_value(75);
  RadiusFrame.set_label("Ion radii");
  RadiusFrame.add(RadiusScale);
  RadiusBox.pack_start(RadiusFrame, Gtk::PACK_SHRINK, 5);
  OptionsBox.pack_start(RadiusBox,  Gtk::PACK_SHRINK, 5);

  ///////////////////////////////////////
  // Setup multiband selection widgets //
  ///////////////////////////////////////
  VisibleBandFrame.add (VisibleBandBox);
  VisibleBandBox.pack_start (MultiBandButton, Gtk::PACK_SHRINK, 5);
  VisibleBandBox.pack_start (VisibleBandWindow);
  VisibleBandWindow.add (VisibleBandTable);
  VisibleBandFrame.set_label("Multiband");
  OptionsBox.pack_start(VisibleBandFrame, Gtk::PACK_SHRINK, 5);
  VisibleBandWindow.property_height_request().set_value(680);
  VisibleBandWindow.set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
  MultiBandButton.set_label ("Enable");
  MultiBandButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnMultiBandToggle));

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
  ((std::vector<Gtk::Widget*>)xPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
  yPlaneButton.set_label("y Plane");
  ((std::vector<Gtk::Widget*>)yPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("green"));
  zPlaneButton.set_label("z Plane");
  ((std::vector<Gtk::Widget*>)zPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("blue"));
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
  Actions->add (Gtk::Action::create("OpenState", "Open State"),
		sigc::mem_fun(*this, &WFVisualClass::OnOpenState));
  Actions->add (Gtk::Action::create("SaveState", "Save State"),
		sigc::mem_fun(*this, &WFVisualClass::OnSaveState));
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

  BoxToggle = Gtk::ToggleAction::create("Box", "Box", "Show box", true);
  Actions->add (BoxToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnBoxToggle));

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


  ///////////////////
  // Colormap menu //
  ///////////////////
  vector<string> mapNames;
  mapNames.push_back ("Autumn");
  mapNames.push_back ("Bone");
  mapNames.push_back ("Colorcube");
  mapNames.push_back ("Cool");
  mapNames.push_back ("Copper");
  mapNames.push_back ("Flag");
  mapNames.push_back ("Gray");
  mapNames.push_back ("Hot");
  mapNames.push_back ("HSV");
  mapNames.push_back ("Jet");
  mapNames.push_back ("Lines");
  mapNames.push_back ("Pink");
  mapNames.push_back ("Spring");
  mapNames.push_back ("Summer");
  mapNames.push_back ("White");
  mapNames.push_back ("Winter");
  mapNames.push_back ("BlueWhiteRed");
  Actions->add (Gtk::Action::create("MenuColormap", "Colormap"));
  CMapActions.resize(BLUE_WHITE_RED+1);
  for (int i=0; i <= BLUE_WHITE_RED; i++) {
    CMapActions[i] = Gtk::RadioAction::create
      (ColorMapGroup, mapNames[i], mapNames[i]);
    Actions->add
      (CMapActions[i], sigc::bind<ColorMapType>
       (sigc::mem_fun (*this, &WFVisualClass::OnColorMapRadio),
	(ColorMapType)i));
  }

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <menuitem action='SaveState'/>"
    "      <menuitem action='OpenState'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
    "      <menuitem action='Nuclei'/>"
    "      <menuitem action='Axes'/>"
    "      <menuitem action='Box'/>"
    "    </menu>"
    "    <menu action='MenuDisplay'>"
    "      <menuitem action='Real'/>"
    "      <menuitem action='Imag'/>"
    "      <menuitem action='Mag2'/>"
    "    </menu>"
    "    <menu action='MenuColormap'>"
    "       <menuitem action='Autumn'/>"
    "       <menuitem action='Bone'/>"
    "       <menuitem action='Colorcube'/>"
    "       <menuitem action='Cool'/>"
    "       <menuitem action='Copper'/>"
    "       <menuitem action='Flag'/>"
    "       <menuitem action='Gray'/>"
    "       <menuitem action='Hot'/>"
    "       <menuitem action='HSV'/>"
    "       <menuitem action='Jet'/>"
    "       <menuitem action='Lines'/>"
    "       <menuitem action='Pink'/>"
    "       <menuitem action='Spring'/>"
    "       <menuitem action='Summer'/>"
    "       <menuitem action='White'/>"
    "       <menuitem action='Winter'/>"
    "       <menuitem action='BlueWhiteRed'/>"
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

  ////////////////////
  // Setup choosers //
  ////////////////////
  SaveStateChooser.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  SaveStateChooser.add_button (Gtk::Stock::SAVE,   Gtk::RESPONSE_OK);
  OpenStateChooser.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  OpenStateChooser.add_button (Gtk::Stock::OPEN,   Gtk::RESPONSE_OK);

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
  MiddleBox.pack_start(PathVis, Gtk::PACK_SHRINK, 5);
  MiddleBox.pack_start(OptionsBox, Gtk::PACK_SHRINK, 5);
  //  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(MiddleBox, Gtk::PACK_SHRINK, 5);
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
  bool boxVisible = BoxToggle->get_active();
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
    Vec3 box = Box(0) + Box(1) + Box(2);
    coord->Set (box);
    PathVis.Objects.push_back(coord);
  }

  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  //boxObject->Set (Box, clipping);
  boxObject->Set (Box.GetLattice(), boxVisible, clipping);
  PathVis.Objects.push_back(boxObject);
  
  Vec3 nLattice[3];
  double length[3];
  length[0] = sqrt(dot(Box(0), Box(0)));
  length[1] = sqrt(dot(Box(1), Box(1)));
  length[2] = sqrt(dot(Box(2), Box(2)));
  nLattice[0] = Box(0)/length[0];
  nLattice[1] = Box(1)/length[1];
  nLattice[2] = Box(2)/length[2];
  
  if (SphereToggle->get_active()) {
    list<AtomClass>::iterator iter;
    list<AtomClass> sphereList;
    for (int ptcl=0; ptcl < AtomPos.extent(0); ptcl++) 
      sphereList.push_back(AtomClass(AtomPos(ptcl), AtomTypes(ptcl)));
    if (clipping) {
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	Vec3 n = Box.GetLatticeInv() * r; 
	int type = (*iter).Type;
	double radius = RadiusScale.get_value() *
	  ElementData::GetRadius(type);
	if ((n[0]+radius/length[0]) > 0.5) 
	  sphereList.push_front(AtomClass(r-Box(0), type));
	if ((n[0]-radius/length[0]) < -0.5) 
	  sphereList.push_front(AtomClass(r+Box(0), type));
      }
      
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	Vec3 n = Box.GetLatticeInv() * r; 
	int type = (*iter).Type;
	double radius = RadiusScale.get_value() *
	  ElementData::GetRadius(type)/length[0];
	if ((n[1]+radius/length[1]) > 0.5)
	  sphereList.push_front(AtomClass(r-Box(1),type));
	if ((n[1]-radius/length[1]) < -0.5)
	  sphereList.push_front(AtomClass(r+Box(1),type));
      }
      
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	Vec3 n = Box.GetLatticeInv() * r; 
	int type = (*iter).Type;
	double radius = RadiusScale.get_value() *
	  ElementData::GetRadius(type);
	if ((n[2]+radius/length[2]) > 0.5)
	  sphereList.push_front(AtomClass(r-Box(2),type));
	if ((n[2]-radius/length[2]) < -0.5)
	  sphereList.push_front(AtomClass(r+Box(2),type));
      }
      // Now make disks to close spheres
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = RadiusScale.get_value() *
	  ElementData::GetRadius(type);
	
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
      double radius = RadiusScale.get_value() *
	ElementData::GetRadius(iter->Type);
      sphere->SetColor(color);
      sphere->SetBox(Box.GetLattice());
      sphere->SetRadius(radius);
      PathVis.Objects.push_back(sphere);
    }
  }

  if (FileIsOpen && !MultiBandButton.get_active()) {
    int band, k;
    band = (int)round(BandScale.get_value());
    k    = (int)round(kScale.get_value());
    if ((CurrBand != band) || (Currk != k)) {
      CurrBand = band;
      Currk    = k;
      ReadWF (k, band);
      Xgrid.Init(-0.5, 0.5, WFData.extent(0));
      Ygrid.Init(-0.5, 0.5, WFData.extent(1));
      Zgrid.Init(-0.5, 0.5, WFData.extent(2));
      WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
      // WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
      WFIso.SetLattice (Box.GetLattice());
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }
    if (ResetIso) {
      Xgrid.Init(-0.5, 0.5, WFData.extent(0));
      Ygrid.Init(-0.5, 0.5, WFData.extent(1));
      Zgrid.Init(-0.5, 0.5, WFData.extent(2));
      WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
      // WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
      WFIso.SetLattice (Box.GetLattice());
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
	vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
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
  if (MultiBandButton.get_active()) {
    if (UpdateIsoVal || UpdateIsoType)
      UpdateMultiIsos();
    for (int i=0; i<VisibleBandRows.size(); i++) {
      BandRow& band = *(VisibleBandRows[i]);
      if (band.Iso != NULL)
	PathVis.Objects.push_back(band.Iso);
    }
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

Mat3 ToMat3 (Array<double,2> &a)
{
  assert (a.rows()==3);
  assert (a.cols()==3);
  Mat3 b;
  b = 
    a(0,0), a(0,1), a(0,2), 
    a(1,0), a(1,1), a(1,2),
    a(2,0), a(2,1), a(2,2);
  return b;
}


void
WFVisualClass::SetupBandTable()
{
  kLabel.set_text ("k");
  BandLabel.set_text ("Band");
  VisibleBandTable.resize(Numk*NumBands+1,3);
  VisibleBandTable.attach (kLabel,    1, 2, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandTable.attach (BandLabel, 2, 3, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandRows.resize (Numk*NumBands);
  int row = 0;
  for (int ki=0; ki<Numk; ki++) 
    for (int bi=0; bi<NumBands; bi++) {
      BandRow &band = *(new BandRow(*this));
      band.Check.signal_toggled().connect 
	(sigc::bind<int>
	 (sigc::mem_fun(*this, &WFVisualClass::OnBandToggle),row));
      VisibleBandTable.attach(band.Check,     0, 1, row+1, row+2, 
			      Gtk::EXPAND, Gtk::SHRINK);
      VisibleBandTable.attach(band.kLabel,    1, 2, row+1, row+2,
			      Gtk::EXPAND, Gtk::SHRINK);
      VisibleBandTable.attach(band.BandLabel, 2, 3, row+1, row+2, 
			      Gtk::EXPAND, Gtk::SHRINK);
      char kstr[50], bstr[50];
      snprintf (kstr, 50, "%d", ki+1);
      snprintf (bstr, 50, "%d", bi+1);
      band.kLabel.set_text(kstr);
      band.BandLabel.set_text(bstr);
      if (bi < NumElectrons/2) {
	band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
	band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
      }
      else {
	band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
	band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
      }
      VisibleBandRows[row] = &band;
      row++;
    }		   
  VisibleBandTable.show_all();
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
  assert (Infile.ReadVar("num_electrons", NumElectrons));
//   if (!IsDiag(lattice)) { 
//     cerr << "wfvis does not currently work with non-orthorhombic cells.\n";
//     abort();
//     //Box.Set (ToMat3(lattice));
//   }
  //  Box.Set (lattice(0,0), lattice(1,1), lattice(2,2));
  Box.Set (ToMat3(lattice));
  xPlane.SetLattice (ToMat3(lattice));
  yPlane.SetLattice (ToMat3(lattice));
  zPlane.SetLattice (ToMat3(lattice));
  
  Infile.CloseSection(); // parameters

  /// Read ion positions
  assert (Infile.OpenSection("ions"));
  Array<double,2> pos;
  assert (Infile.ReadVar("pos", pos));
  AtomPos.resize(pos.extent(0));
//   for (int i=0; i<pos.extent(0); i++)
//     for (int j=0; j<3; j++)
//       AtomPos(i)[j] = (pos(i,j) - 0.5*Box[j]);
  for (int i=0; i<pos.extent(0); i++) {
    AtomPos(i) = Vec3(pos(i,0), pos(i,1), pos(i,2));
    for (int j=0; j<3; j++)
      AtomPos(i) -= (0.5-Shift[j])*Box(j);
  }
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
  SetupBandTable();

  /// Read first wave function
  ReadWF(0,0);
//   Xgrid.Init(-0.5*Box[0], 0.5*Box[0], WFData.extent(0));
//   Ygrid.Init(-0.5*Box[1], 0.5*Box[1], WFData.extent(1));
//   Zgrid.Init(-0.5*Box[2], 0.5*Box[2], WFData.extent(2));
  Xgrid.Init(-0.5, 0.5, WFData.extent(0));
  Ygrid.Init(-0.5, 0.5, WFData.extent(1));
  Zgrid.Init(-0.5, 0.5, WFData.extent(2));
  WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
  //  WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
  WFIso.SetLattice(Box.GetLattice());
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
  UpdateIsoVal = true;
  DrawFrame();
}

void
WFVisualClass::OnBandChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
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
  Array<double,1> twist_angle;
  assert (Infile.ReadVar("twist_angle", twist_angle));
  bool gammaPoint = 
    (fabs(twist_angle(0)) < 1.0e-12) &&
    (fabs(twist_angle(1)) < 1.0e-12) &&
    (fabs(twist_angle(2)) < 1.0e-12);
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
  // If we're at the gamma point, adjust the phase of the WF so that
  // it is purely real.
  if (gammaPoint && ((WFDisplay == REAL_PART) || (WFDisplay == IMAG_PART))) {
    double maxRho = 0.0;
    int ixMax, iyMax, izMax;
    for (int ix=0; ix<Nx; ix+=5)
      for (int iy=0; iy<Ny; iy+=5)
	for (int iz=0; iz<Nz; iz+=5) {
	  double rho = (wfdata(ix,iy,iz,0)*wfdata(ix,iy,iz,0) +
			wfdata(ix,iy,iz,1)*wfdata(ix,iy,iz,1));
	  if (rho > maxRho) {
	    maxRho = rho;
	    ixMax = ix; iyMax = iy; izMax=iz;
	  }
	}
    double phase = atan2 (wfdata(ixMax, iyMax, izMax,1),
			  wfdata(ixMax, iyMax, izMax,0));
    complex<double> factor (cos(phase), -sin(phase));
    for (int ix=0; ix<wfdata.extent(0); ix++)
      for (int iy=0; iy<wfdata.extent(1); iy++)
	for (int iz=0; iz<wfdata.extent(2); iz++) {
	  complex<double> wf(wfdata(ix,iy,iz,0), wfdata(ix,iy,iz,1));	  
	  wf *= factor;
	  wfdata(ix,iy,iz,0) = wf.real();
	  wfdata(ix,iy,iz,1) = wf.imag();
	}
  }
    
  MaxVal = 0.0;
  int xShift, yShift, zShift;
  xShift = (int)round(Shift[0]*wfdata.extent(0));
  yShift = (int)round(Shift[1]*wfdata.extent(1));
  zShift = (int)round(Shift[2]*wfdata.extent(2));
  for (int ix=0; ix<wfdata.extent(0); ix++)
    for (int iy=0; iy<wfdata.extent(1); iy++)
      for (int iz=0; iz<wfdata.extent(2); iz++) {
	// WF data is store from 0 to Lx, not -Lx/2 to Lx/2
	int jx = (ix-xShift+Nx-1/*+Nx/2*/)%(Nx-1);
	int jy = (iy-yShift+Ny-1/*+Ny/2*/)%(Ny-1);
	int jz = (iz-zShift+Nz-1/*+Nz/2*/)%(Nz-1);
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
  cerr << "MaxVal = " << MaxVal << endl;
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
WFVisualClass::OnBoxToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnRadiusChange()
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
  UpdateIsoType = true;
  if (MultiBandButton.get_active())
    DrawFrame();
}


bool
WFVisualClass::WriteState(string fname)
{
  IOSectionClass out;
  
  if (out.NewFile(fname) == false)
    return false;

  out.NewSection ("Flags");
  out.WriteVar ("xPlane", xPlaneButton.get_active());
  out.WriteVar ("yPlane", yPlaneButton.get_active());
  out.WriteVar ("zPlane", zPlaneButton.get_active());
  out.WriteVar ("Nuclei", SphereToggle->get_active());
  out.WriteVar ("CoordAxes", CoordToggle->get_active());
  out.WriteVar ("Isosurface", IsoButton.get_active());
  out.WriteVar ("Clip", ClipButton.get_active());
  out.WriteVar ("Perspective", PerspectButton.get_active());
  out.WriteVar ("MultiBands", MultiBandButton.get_active());
  Array<bool,1> visibleBands (VisibleBandRows.size());
  for (int i=0; i<VisibleBandRows.size(); i++)
    visibleBands(i) = VisibleBandRows[i]->Check.get_active();
  out.WriteVar ("VisibleBands", visibleBands);

  string wftype;
  if (WFDisplay == REAL_PART)
    wftype = "real";
  if (WFDisplay == IMAG_PART)
    wftype = "imag";
  if (WFDisplay == MAG2)
    wftype = "mag2";
  out.WriteVar ("WFDisplay", wftype);
  out.CloseSection();

  out.WriteVar ("kPoint", (int)round(kAdjust.get_value()));
  out.WriteVar ("Band", (int)round(BandAdjust.get_value()));
  out.WriteVar ("IsoVal", IsoAdjust.get_value());
  out.WriteVar ("xPlanePos", xPlaneAdjust.get_value());
  out.WriteVar ("yPlanePos", yPlaneAdjust.get_value());
  out.WriteVar ("zPlanePos", zPlaneAdjust.get_value());
  out.WriteVar ("Radius", RadiusScale.get_value());
  out.NewSection("View");
  PathVis.View.Write(out);
  out.CloseSection();
  out.CloseFile();
  return true;
}

bool
WFVisualClass::ReadState (string fname)
{
  IOSectionClass in;
  if (in.OpenFile(fname) == false)
    return false;
  bool active;
  
  assert(in.OpenSection ("Flags"));
  in.ReadVar ("xPlane", active);
  xPlaneButton.set_active(active);
  assert(in.ReadVar ("yPlane", active));
  yPlaneButton.set_active(active);
  assert(in.ReadVar ("zPlane", active));
  zPlaneButton.set_active(active);
  assert(in.ReadVar ("Nuclei", active));
  SphereToggle->set_active(active);
  assert(in.ReadVar ("Isosurface", active));
  IsoButton.set_active(active);
  assert(in.ReadVar ("Clip", active));
  ClipButton.set_active(active);
  assert(in.ReadVar ("Perspective", active));
  PerspectButton.set_active(active);
  if (in.ReadVar ("CoordAxes", active))
    CoordToggle->set_active(active);
  if (in.ReadVar ("MultiBands", active))
    MultiBandButton.set_active(active);
  Array<bool,1> visibleBands;
  if (in.ReadVar ("VisibleBands", visibleBands))
    for (int i=0; i<visibleBands.size(); i++) 
      VisibleBandRows[i]->Check.set_active(visibleBands(i));


  string wftype;
  assert(in.ReadVar ("WFDisplay", wftype));
  if (wftype == "real") {
    RealRadio->set_active(true);
    WFDisplay = REAL_PART;
  }
  else if (wftype == "imag") {
    ImagRadio->set_active(true);
    WFDisplay = IMAG_PART;
  }
  else if (wftype == "mag2") {
    Mag2Radio->set_active(true);
    WFDisplay = MAG2;
  }
  in.CloseSection(); // flags

  int k, band;
  assert(in.ReadVar ("kPoint", k));
  assert(in.ReadVar ("Band", band));
  double val;
  assert(in.ReadVar ("IsoVal", val));
  IsoAdjust.set_value(val);

  assert(in.ReadVar ("xPlanePos", val));
  xPlaneAdjust.set_value(val);
  assert(in.ReadVar ("yPlanePos", val));
  yPlaneAdjust.set_value(val);
  assert(in.ReadVar ("zPlanePos", val));
  zPlaneAdjust.set_value(val);
  if (in.ReadVar("Radius", val))
    RadiusScale.set_value(val);

  kAdjust.set_value((double)k);
  BandAdjust.set_value((double)band);
  assert(in.OpenSection("View"));
  PathVis.View.Read(in);
  in.CloseSection();
  in.CloseFile();
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
  return true;
}

void
WFVisualClass::OnSaveState()
{
  int result = SaveStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = SaveStateChooser.get_filename();
    WriteState (fname);
  }
  SaveStateChooser.hide();
}

void
WFVisualClass::OnOpenState()
{
  int result = OpenStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = OpenStateChooser.get_filename();
    ReadState (fname);
  }
  OpenStateChooser.hide();
}  



void
WFVisualClass::OnMultiBandToggle()
{
  if (MultiBandButton.get_active()) {
    IsoButton.set_active(false);
    xPlaneButton.set_active(false);
    yPlaneButton.set_active(false);
    zPlaneButton.set_active(false);
    IsoButton.set_sensitive (false);
    PlaneFrame.set_sensitive(false);
    BandFrame.set_sensitive(false);
    kFrame.set_sensitive(false);
  }
  else {
    IsoButton.set_sensitive (true);
    PlaneFrame.set_sensitive(true);
    BandFrame.set_sensitive (true);
    kFrame.set_sensitive    (true);
  }
}

void
WFVisualClass::OnBandToggle (int row)
{
  int ki = row/NumBands;
  int bi = row % NumBands;
  BandRow &band = *(VisibleBandRows[row]);
  if (band.Check.get_active()) {
    if (band.Iso == NULL) {
      band.Iso = new Isosurface;
      if (FileIsOpen) {
	ReadWF (ki, bi);
	Xgrid.Init(-0.5, 0.5, WFData.extent(0));
	Ygrid.Init(-0.5, 0.5, WFData.extent(1));
	Zgrid.Init(-0.5, 0.5, WFData.extent(2));
	band.Iso->Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
	//	band.Iso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
	band.Iso->SetLattice (Box.GetLattice());
	if (WFDisplay == MAG2) {
	  band.Iso->SetColor (0.0, 0.8, 0.0);
	  band.Iso->SetIsoval(MaxVal*IsoAdjust.get_value());
	}
	else {
	  vector<TinyVector<double,3> > colors;
	  colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	  colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	  band.Iso->SetColor(colors);
	  vector<double> vals;
	  vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	  vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
	  band.Iso->SetIsoval(vals);
	}
      }
      band.Iso->Dynamic = false;
    } 
  }
  else {
    if (band.Iso != NULL) {
      delete band.Iso;
      band.Iso = NULL;
    }
  }
  if (MultiBandButton.get_active())
    DrawFrame();
}

void
WFVisualClass::UpdateMultiIsos()
{
  for (int i=0; i<VisibleBandRows.size(); i++) {
    int ki = i / NumBands;
    int bi = i % NumBands;
    if (VisibleBandRows[i]->Iso != NULL) {
      Isosurface &iso = *(VisibleBandRows[i]->Iso);
      if (UpdateIsoType) 
	ReadWF(ki, bi);
      if (UpdateIsoType || UpdateIsoVal) {
	if (WFDisplay == MAG2) {
	  iso.SetColor (0.0, 0.8, 0.0);
	  iso.SetIsoval(MaxVal*IsoAdjust.get_value());
	}
	else {
	  vector<TinyVector<double,3> > colors;
	  colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	  colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	  iso.SetColor(colors);
	  vector<double> vals;
	  vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	  vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
	  iso.SetIsoval(vals);
	}
      }
    }
  }
  UpdateIsoVal = false;
  UpdateIsoType = false;
}

void
WFVisualClass::SetShift(Vec3 shift)
{
  DoShift = true;
  Shift = shift;
}


void
WFVisualClass::SetViewportSize (int size)
{
  VisibleBandWindow.property_height_request().set_value(size-100);
  PathVis.set_size_request(size, size);
  resize(10,10);
}

void
WFVisualClass::OnColorMapRadio (ColorMapType type)
{
  if (CMapActions[type]->get_active()) {
    CMap = type;
    xPlane.SetColorMap(type);
    yPlane.SetColorMap(type);
    zPlane.SetColorMap(type);
    DrawFrame();
  }
}


vector<string>
BreakString (string str, char sep)
{
  vector<string> strvec;
  int len = str.size();
  int i = 0;
  while (i < len) {
    char s[2];
    s[1] = '\0';
    string item;
    while (str[i] != sep && i < len) {
      s[0] = str[i];
      item.append(s);
      i++;
    }
    strvec.push_back(item);
    i++;
  }
  return strvec;
}


int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);
  glutInit(&argc, argv);

  list<ParamClass> optionList;
  optionList.push_back(ParamClass("shift", true));
  optionList.push_back(ParamClass("small", false));
  optionList.push_back(ParamClass("remote", false));
  CommandLineParserClass parser (optionList);
  bool success = parser.Parse (argc, argv);
  if (!success || parser.NumFiles() < 1 || parser.NumFiles() > 2) {
    cerr << "Usage:\n  wfvis++ [options...] myfile.h5 [statefile.h5]\n"
	 << "Options:\n"
	 << "  --shift x,y,z       shift by reduced coordinates\n"
	 << "  --small             reduce size for small displays\n"
	 << "  --remote            reduce data transfer for remote operation\n"
	 << "                      over a network connection\n";
    exit (1);
  }
  

  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  // Instantiate and run the application.
  WFVisualClass wfvisual;

  if (parser.Found("shift")) {
    string shiftStr = parser.GetArg("shift");
    Vec3 shift;
    vector<string> components = BreakString (shiftStr, ',');
    if (components.size() != 3) {
      cerr << "Expected 3 components for shift.\n";
      abort();
    }
    for (int i=0; i<3; i++)
      shift[i] = strtod (components[i].c_str(), NULL);
    wfvisual.SetShift (shift);
  }
  if (parser.Found("small"))
    wfvisual.SetViewportSize (600);

  wfvisual.Read (parser.GetFile(0));

  if (parser.NumFiles() == 2)
    wfvisual.ReadState (parser.GetFile(1));
  kit.run(wfvisual);

  return 0;
}
