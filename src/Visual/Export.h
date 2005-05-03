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
  /// These are used for GL-based exporting
  Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap;
  Glib::RefPtr<Gdk::Pixmap>     GdkPixmap;
  Glib::RefPtr<Gdk::Drawable>   Draw;
  Glib::RefPtr<Gdk::GL::Config> GLConfig;
  Glib::RefPtr<Gdk::GL::Context> GLContext;
  /// Widgets for the UI
  Gtk::SpinButton WidthButton, HeightButton;
  Gtk::HBox WidthBox, HeightBox;
  Gtk::Label WidthLabel, HeightLabel;
  Gtk::Frame TypeFrame;
  Gtk::HBox TypeBox;
  Gtk::ComboBoxText TypeCombo;
  Gtk::HBox WidthHeightBox;
  Gtk::VBox MainVBox;
  Gtk::HButtonBox ButtonBox;
  Gtk::Button ExportButton, CancelButton, BrowseButton;
  Gtk::Frame BaseNameFrame;
  Gtk::HBox BaseNameHBox;
  Gtk::Entry BaseNameEntry;
  Gtk::Button BaseNameBrowseButton;
  Gtk::FileChooserDialog BaseNameChooser;
  /// Signal handlers
  void OnExportButton(), OnCancelButton(), OnBrowseButton();

  int Width, Height;
  int AttribList;
  void InitGLStuff();
public:
  void Export (string filename);
  void ExportPOV (string filename);
  ExportClass (VisualClass &visual) : 
    Visual(visual), Width(2000), Height(2000),
    BaseNameChooser("Export filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
    WidthLabel("Width:"), HeightLabel("Height:")
  {
    BaseNameFrame.set_label("Base filename");
    BaseNameFrame.add (BaseNameHBox);
    BaseNameHBox.pack_start(BaseNameEntry, Gtk::PACK_SHRINK, 5);
    BaseNameHBox.pack_start(BaseNameBrowseButton, Gtk::PACK_SHRINK);
    BaseNameBrowseButton.set_label("Browse");
    MainVBox.pack_start (BaseNameFrame, Gtk::PACK_SHRINK, 5);
    WidthBox.pack_start (WidthLabel, Gtk::PACK_SHRINK,5);
    WidthBox.pack_start (WidthButton, Gtk::PACK_SHRINK);
    HeightBox.pack_start (HeightLabel, Gtk::PACK_SHRINK,5);
    HeightBox.pack_start (HeightButton, Gtk::PACK_SHRINK);
    WidthHeightBox.pack_start (WidthBox, Gtk::PACK_SHRINK, 5);
    WidthHeightBox.pack_start (HeightBox, Gtk::PACK_SHRINK, 5);
    TypeFrame.set_label("Rendering Type");
    TypeFrame.add(TypeBox);
    TypeBox.pack_start(TypeCombo, Gtk::PACK_SHRINK);
    TypeCombo.insert_text(0, "Offscreen GL");
    TypeCombo.insert_text(1, "POVray");
    TypeCombo.set_active (0);
    ExportButton.set_label ("Export");
    CancelButton.set_label ("Cancel");
    ButtonBox.pack_start(ExportButton,  Gtk::PACK_SHRINK, 5);
    ButtonBox.pack_start(CancelButton,  Gtk::PACK_SHRINK, 5);
    MainVBox.pack_start(WidthHeightBox, Gtk::PACK_SHRINK, 5);
    MainVBox.pack_start(TypeFrame,      Gtk::PACK_SHRINK, 5);
    MainVBox.pack_start(ButtonBox, Gtk::PACK_SHRINK, 10);
    add (MainVBox);
    set_title ("Export");
    // Set signal handlers
    ExportButton.signal_clicked().connect
      (sigc::mem_fun(*this, &ExportClass::OnExportButton));
    CancelButton.signal_clicked().connect
      (sigc::mem_fun(*this, &ExportClass::OnCancelButton));
    BaseNameBrowseButton.signal_clicked().connect
      (sigc::mem_fun(*this, &ExportClass::OnBrowseButton));
   };
};


#endif
