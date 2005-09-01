#!/usr/bin/python 
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.numerix import arange, sin, pi
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import gtk
import math
from IO import *

class HelloWorld:
    def delete_event(self, widget, event, data=None):
        # Do any popup checks etc.
        return False

    def destroy(self, widget, data=None):
        gtk.main_quit()

    def create_tree(self, name):
        self.treestore = gtk.TreeStore(str)
#        for parent in range(4):
#            piter = self.treestore.append(None, ['parent %i' % parent])
#            for child in range(3):
#                self.treestore.append(piter,['child %i of parent %i' % (child, parent)])

        self.treeview = gtk.TreeView(self.treestore)
        self.tvcolumn = gtk.TreeViewColumn('Observables')
        self.treeview.append_column(self.tvcolumn)
        self.cell = gtk.CellRendererText()
        self.tvcolumn.pack_start(self.cell, True)
        self.tvcolumn.add_attribute(self.cell, 'text', 0)
        self.treeview.set_search_column(0)
        self.treeview.set_reorderable(True)
        self.vbox.pack_start(self.treeview)
        


    def __init__(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.destroy)

        self.window.set_border_width(10)

        self.vbox = gtk.VBox()
        self.window.add(self.vbox)

        # Create an ActionGroup named ActionExample
        actiongroup = gtk.ActionGroup('MyActionGroup') 

        # Create an action for quitting the program using a stock item
        quitaction = gtk.Action('Quit', '_Quit me!', 'Quit the Program',
                            gtk.STOCK_QUIT)
        quitaction.set_property('short-label', '_Quit')
        # Connect a callback to the action
        quitaction.connect('activate', self.quit_cb)

        # Create an action for openting the program using a stock item
        openaction = gtk.Action('Open', '_Open me!', 'Open the Program',
                            gtk.STOCK_OPEN)
        openaction.set_property('short-label', '_Open')
        # Connect a callback to the action
        openaction.connect('activate', self.open_cb)



        # Create the MenuBar
        menubar = gtk.MenuBar()
        self.vbox.pack_start(menubar, False)

        # Create File menu stuff
        file_action = gtk.Action('File', '_File', None, None)
        actiongroup.add_action(file_action)
        file_menuitem = file_action.create_menu_item()
        menubar.append(file_menuitem)
        file_menu = gtk.Menu()
        file_menuitem.set_submenu(file_menu)

        # Create a proxy MenuItem
        openitem = openaction.create_menu_item()
        quititem = quitaction.create_menu_item()
        file_menu.append(openitem)
        file_menu.append(quititem)
        
        # Create a Toolbar
        toolbar = gtk.Toolbar()
        self.vbox.pack_start(toolbar, False)

        # Create a proxy ToolItem
        opentoolitem = openaction.create_tool_item()
        quittoolitem = quitaction.create_tool_item()
        toolbar.insert(opentoolitem, 0)
        toolbar.insert(quittoolitem, 0)
        self.create_tree('abc')
        self.window.set_size_request(400,500)
        self.window.show_all()

        self.chooser = gtk.FileChooserDialog(title='Select a file', action=gtk.FILE_CHOOSER_ACTION_OPEN, buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        filter = gtk.FileFilter()
        filter.set_name ("HDF5 files")
        filter.add_pattern("*.h5")
        self.chooser.add_filter(filter)
        filter = gtk.FileFilter()
        filter.set_name ("All files")
        filter.add_pattern("*")
        self.chooser.add_filter(filter)

    def quit_cb(self, b):
        print 'Quitting program'
        gtk.main_quit()

    def open_file (self, name):
        print 'Opening file ' + name
        self.infile = IOSectionClass()
        success = self.infile.OpenFile(name)
        if (success != True):
            print 'Cannot open file ' + name
        else:
            print 'Opened file successfully'
            self.infile.OpenSection("Observables")
            numSecs = self.infile.CountSections()
            gofrRows = self.treestore.append(None, ["g(r)"])
            for i in range(numSecs):
                self.infile.OpenSection(i)
                myType = self.infile.ReadVar("Type")

                if (myType == "CorrelationFunction"):
                    iter = self.treestore.append(gofrRows, [self.infile.GetName()])
                    gofr = gofrClass(self.infile, self.treestore.get_path(iter))
                    self.treeview.connect_object("row-activated", gofrClass.tv_callback, gofr)
                elif (myType == "Scalar"):
                    iter = self.treestore.append(None, [self.infile.GetName()])
                    numVars = self.infile.CountVars()
                    for j in range (numVars):
                        varName = self.infile.GetVarName(j)
                        if ((varName != "Description") & (varName != "Type")):
                            variter = self.treestore.append (iter, [varName])
                            scalar = scalarClass(self.infile, varName, self.treestore.get_path(variter))
                            self.treeview.connect_object("row-activated", scalarClass.tv_callback, scalar)
                elif (myType == "Path"):
                    iter = self.treestore.append(None, ["Path Dump"])
                self.infile.CloseSection()
            self.infile.CloseSection() # Observables
        return True
        
    def open_cb(self, b):
        response = self.chooser.run()
        if response == gtk.RESPONSE_OK:
            self.open_file (self.chooser.get_filename())
        self.chooser.hide()
        return False


class ScalarVar:
    def __init__(self, infile, varname):
        self.name = varname

class gofrClass:
    def __init__(self, infile, rowNum):
        self.rowNum = rowNum
        self.name = infile.GetName()
        self.x = infile.ReadVar("x")
        self.y = infile.ReadVar("y")
        if (infile.ReadVar("Cumulative") == False):
            n = self.y.shape[0]
            self.y = sum(self.y) / n
        print "New g(r) class with name " + self.name

    def tv_callback (self, path, view_column):
        if (path == self.rowNum):
            self.plot()
    
    def plot(self):
        win = gtk.Window()
        vbox = gtk.VBox()
        win.add (vbox)
        fig = Figure(figsize=(5,4), dpi=100)
        ax = fig.add_subplot(111)
        ax.plot(self.x,self.y)
        ax.set_title(self.name)
        ax.set_xlabel("r")
        ax.set_ylabel(r'$g(r)$')
#        ax.title = self.name
        canvas = FigureCanvas(fig)
        vbox.pack_start(canvas)
        toolbar = NavigationToolbar(canvas, win)
        vbox.pack_start(toolbar, False, False)
        win.set_size_request(520,500)
        win.show_all()

class scalarClass:
    def __init__(self, infile, varName, rowNum):
        self.rowNum = rowNum
        self.name = varName
        self.y = infile.ReadVar(varName)
        self.x = arange(0,self.y.size())
        print "New scalar class with name " + self.name

    def tv_callback (self, path, view_column):
        if (path == self.rowNum):
            self.plot()
    
    def plot(self):
        win = gtk.Window()
        vbox = gtk.VBox()
        win.add (vbox)
        fig = Figure(figsize=(5,4), dpi=100)
        ax = fig.add_subplot(111)
        ax.plot(self.x,self.y)
        ax.set_title(self.name)
        ax.set_xlabel("block")
        ax.set_ylabel("Energy")
        canvas = FigureCanvas(fig)
        vbox.pack_start(canvas)
        toolbar = NavigationToolbar(canvas, win)
        vbox.pack_start(toolbar, False, False)
        win.set_size_request(520,500)
        win.show_all()

if __name__ == '__main__':
    ba = HelloWorld()
    gtk.main()



