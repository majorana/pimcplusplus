from distutils.core import setup, Extension

module1 = Extension('IOSection',
                    sources = ['InputOutputPythonWrapper.cc','InputOutput.cc',\
                               'InputOutputHDF5.cc', 'InputOutputASCII.cc',\
                               'InputOutputXML.cc' ],\
                    include_dirs=['/usr/include/libxml2'],
                    library_dirs=[ ],
                    libraries =  ['blitz', 'xml2', 'hdf5'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])
