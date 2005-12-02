#include "MDVis.h"

MDVisualClass::MDVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  FrameAdjust(0.0, 0.0, 0.0),
  CurrentFrame(0),
  PlayDirection(1),
  TimeoutDelay(20)
{
  FrameScale.set_adjustment(FrameAdjust);
  FrameScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnFrameChange));

  PlayImage.property_file().set_value(FindFullPath("player_play.png"));
  PlayButton.set_icon_widget(PlayImage);
  PlayButton.set_label("Play");
  RevImage.property_file().set_value (FindFullPath("player_rev.png"));
  RevButton.set_icon_widget(RevImage);
  RevButton.set_label("Reverse");
  PauseImage.property_file().set_value (FindFullPath("player_pause.png"));
  PauseButton.set_icon_widget(PauseImage);
  PauseButton.set_label("Pause");

  OrthoImage.property_file().set_value(FindFullPath("orthographic.png"));
  OrthoButton.set_icon_widget(OrthoImage);
  OrthoButton.set_label("Ortho");
  PerspectImage.property_file().set_value(FindFullPath("perspective.png"));
  PerspectButton.set_icon_widget(PerspectImage);
  PerspectButton.set_label("Perspect");

  set_reallocate_redraws(true);
  PathVis.set_size_request(800, 800);
  ////////////////////
  // Setup tool bar //
  ////////////////////
  Gtk::RadioButtonGroup group =PauseButton.get_group();
  RevButton.set_group(group);
  group = RevButton.get_group();
  PlayButton.set_group(group);

  group = OrthoButton.get_group();
  PerspectButton.set_group(group);
  
  Tools.append(RevButton);
  Tools.append(PauseButton);
  Tools.append(PlayButton);
  Tools.append(ToolSep1);
  Tools.append(OrthoButton);
  Tools.append(PerspectButton);
  
  /////////////////
  // Setup menus //
  /////////////////
  Actions = Gtk::ActionGroup::create();
  Actions->add (Gtk::Action::create("MenuFile", "_File"));
  Actions->add (Gtk::Action::create("Open", "_Open"),
		sigc::mem_fun(*this, &MDVisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &MDVisualClass::OnExport));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &MDVisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &MDVisualClass::OnViewReset));

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
  PlayButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  RevButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  PauseButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  OrthoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPerspectiveToggle));
  QuitButton.signal_clicked().connect
    (sigc::mem_fun(*this, &MDVisualClass::Quit));

  ////////////////////
  // Pack the boxes //
  ////////////////////
  ToolBox.pack_start(Tools);
  MainVBox.pack_start(*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK,0);
  MainVBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(FrameScale, Gtk::PACK_SHRINK, 0);
  MainVBox.pack_start(QuitButton, Gtk::PACK_SHRINK, 0);

  add (MainVBox);
  set_title ("mdvis++");
  show_all();
}

void
MDVisualClass::OnViewReset()
{

}


void
MDVisualClass::OnOpen()
{

}

void
MDVisualClass::Quit()
{
  Gtk::Main::quit();
}

void
MDVisualClass::OnExport()
{

}





string 
MDVisualClass::FindFullPath(string filename)
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
MDVisualClass::OnPlayToggle()
{
  TimeoutConnection.disconnect();
  if (PlayButton.get_active()) 
    PlayDirection = 1;
  else if (RevButton.get_active()) {
    PlayDirection = -1;
  }
  else
    return;
  
  /// Connect idle timeout
  sigc::slot<bool> slot = sigc::mem_fun(*this, &MDVisualClass::OnTimeout);
  TimeoutConnection = Glib::signal_timeout().connect (slot, TimeoutDelay);

}

bool
MDVisualClass::OnTimeout()
{
  DrawFrame();
  PathVis.Invalidate();
  if (PlayDirection == 1) {
    if (CurrentFrame < (Trajectory.extent(0)-1)) {
      CurrentFrame++;
      FrameAdjust.set_value(CurrentFrame);
      return true;
    }
    else {
      PauseButton.set_active();
      return false;
    }
  }
  else if (PlayDirection == -1) {
    if (CurrentFrame > 0) {
      CurrentFrame--;
      FrameAdjust.set_value(CurrentFrame);
      return true;
    }
    else {
      PauseButton.set_active();
      return false;
    }
  }
}

void
MDVisualClass::OnPerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  cerr << "Now using " << (persp ? "perspective" : "orthographic") 
       << " projection.\n";
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  DrawFrame();
  PathVis.Invalidate();
}

void
MDVisualClass::OnFrameChange()
{
  CurrentFrame = (int)round(FrameScale.get_value());
  DrawFrame();
  PathVis.Invalidate();
}

void
MDVisualClass::DrawFrame()
{
  for (int i=0; i<PathVis.Objects.size(); i++)
    delete PathVis.Objects[i];
  PathVis.Objects.resize(0);
  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  boxObject->Set (Box[0], Box[1], Box[2]);
  PathVis.Objects.push_back(boxObject);
  for (int ptcl=0; ptcl<Trajectory.extent(1); ptcl++) {
    SphereObject *sphere = new SphereObject;
    sphere->SetPos(Vec3(Trajectory(CurrentFrame,ptcl,0),
			Trajectory(CurrentFrame,ptcl,1),
			Trajectory(CurrentFrame,ptcl,2)));
    sphere->SetRadius(0.5);
    PathVis.Objects.push_back(sphere);
  }
			
}



void
MDVisualClass::Read(string filename)
{
  IOSectionClass in;
  assert (in.OpenFile(filename));
  Array<double,1> box;
  assert (in.OpenSection("System"));
  assert (in.ReadVar ("Box", box));
  Box.Set (box(0), box(1), box(2));
  in.CloseSection();

  double maxDim = max(max(box(0), box(1)), box(2));
  PathVis.View.SetDistance (1.2*maxDim);
  
  assert(in.OpenSection("Moves"));
  assert(in.OpenSection("Langevin"));
  assert(in.OpenSection("Langevin"));
  assert(in.ReadVar("R", Trajectory));
  assert(in.ReadVar("Time", Time));

  in.CloseSection(); // "Langevin"
  in.CloseSection(); // "Langevin"
  in.CloseSection(); // "Moves"

  in.CloseFile();

  FrameAdjust.set_upper(Trajectory.extent(0)-1);
  FrameScale.set_digits(0);

  DrawFrame();
}


MDVisualClass::~MDVisualClass()
{

}

int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);

  if (argc < 2) {
    cerr << "Usage:\n  mdvis++ myfile.h5\n";
    exit (1);
  }
  
  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  // Instantiate and run the application.
  MDVisualClass mdvisual;
  mdvisual.Read (argv[1]);
  kit.run(mdvisual);

  return 0;
}
