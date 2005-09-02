#!/usr/bin/python 
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.numerix import arange, sin, pi, array
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import sys
import gtk
import math
from IO import *
import stats

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
#        self.tvcolumn.add_attribute(self.cell, 'text', 0)
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
        quitaction = gtk.Action('Quit', '_Quit', 'Quit the Program',
                            gtk.STOCK_QUIT)
        quitaction.set_property('short-label', '_Quit')
        # Connect a callback to the action
        quitaction.connect('activate', self.quit_cb)

        # Create an action for openting the program using a stock item
        openaction = gtk.Action('Open', '_Open', 'Open the Program',
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
        gtk.main_quit()

    def open_file (self, name):
#        print 'Opening file ' + name
        self.infile = IOSectionClass()
        success = self.infile.OpenFile(name)
        if (success != True):
            print 'Cannot open file ' + name
        else:
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
    def __init__(self, infile, path):
        self.path = path
        self.name = infile.GetName()
        self.x = infile.ReadVar("x")
        self.y = infile.ReadVar("y")
        if (infile.ReadVar("Cumulative") == False):
            n = self.y.shape[0]
            self.y = sum(self.y) / n
#        print "New g(r) class with name " + self.name

    def tv_callback (self, path, view_column):
        if (path == self.path):
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
        (self.mean, self.var, self.error, self.kappa) = stats.Stats(self.y)
#        print "New scalar class with name " + self.name

    def tv_callback (self, path, view_column):
        if (path == self.rowNum):
            self.plot()
    
    def plot(self):
        win = gtk.Window()
        winBox = gtk.HBox()
        plotBox = gtk.VBox()
        rightBox = gtk.VBox()
        statFrame = gtk.Frame()
        statBox = gtk.VBox()
        winBox.pack_start(plotBox)
        rightBox.pack_start(statFrame, True, False)
        winBox.pack_start(rightBox, False, False)
#        statFrame.add (statBox)
        win.add (winBox)
        fig = Figure(figsize=(5,4), dpi=100)
        ax = fig.add_subplot(111)
        ax.plot(self.x,self.y)
        ax.hold(True)
        xMean = array([0.0, self.y.size()]) 
        yMean = array([self.mean, self.mean])
        ax.plot(xMean, yMean, 'g')
        ax.set_title(self.name)
        ax.set_xlabel("block")
        ax.set_ylabel("Energy")
        canvas = FigureCanvas(fig)
        plotBox.pack_start(canvas)
        toolbar = NavigationToolbar(canvas, win)
        plotBox.pack_start(toolbar, False, False)
        # Make statistics box
        statFrame.set_label("Statistics")
        statTable = gtk.Table(3, 2, False)
        statFrame.add(statTable)

        (meanStr, errorStr) = MeanErrorString(self.mean, self.error)
        meanLabel1 = gtk.Label("Mean: ")
        meanLabel2 = gtk.Label(meanStr + " +/- " + errorStr)
        statTable.attach(meanLabel1, 0, 1, 0, 1, gtk.SHRINK)
        statTable.attach(meanLabel2, 1, 2, 0, 1, gtk.SHRINK)

        varStr   = '%1.2f' % self.var
        varLabel1 = gtk.Label("Variance:")
        varLabel2 = gtk.Label(varStr)
        statTable.attach(varLabel1, 0, 1, 1, 2, gtk.SHRINK)
        statTable.attach(varLabel2, 1, 2, 1, 2, gtk.SHRINK)

        kappaStr= '%1.2f' % self.kappa
        kappaLabel1=gtk.Label("Kappa: ")
        kappaLabel2=gtk.Label(kappaStr)
        statTable.attach(kappaLabel1, 0, 1, 2, 3, gtk.SHRINK)
        statTable.attach(kappaLabel2, 1, 2, 2, 3, gtk.SHRINK)
        
        win.set_size_request(650,500)
        win.show_all()

def MeanErrorString (mean, error):
     if (mean!=0.0):
          meanDigits = math.floor(math.log(abs(mean))/math.log(10))
     else:
          meanDigits=2
     if (error!=0.0):
          rightDigits = -math.floor(math.log(error)/math.log(10))+1
     else:
          rightDigits=2
     if (rightDigits < 0):
          rightDigits = 0
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return (meanstr, errorstr)

if __name__ == '__main__':
    fileName = sys.argv[1]
    ba = HelloWorld()
    if fileName!="":
        ba.open_file(fileName)
    gtk.main()
    


