from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.numerix import arange, sin, pi

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
        for parent in range(4):
            piter = self.treestore.append(None, ['parent %i' % parent])
            for child in range(3):
                self.treestore.append(piter,['child %i of parent %i' % (child, parent)])

        self.treeview = gtk.TreeView(self.treestore)
        self.tvcolumn = gtk.TreeViewColumn('Column 0')
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
        self.window.show_all()


    def quit_cb(self, b):
        print 'Quitting program'
        gtk.main_quit()

    def open_cb(self, b):
        print 'Open callback'
        return True


if __name__ == '__main__':
    ba = HelloWorld()
    gtk.main()



