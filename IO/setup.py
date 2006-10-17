from distutils.core import setup, Extension

module1 = Extension('IOSection',
                    sources = ['InputOutputPythonWrapper.cc','InputOutput.cc',\
                               'InputOutputHDF5.cc', 'InputOutputASCII.cc',\
                               'InputOutputXML.cc' ],\
                    include_dirs=['/usr/include/libxml2',\
                                  '/home/kesler/include'],
                    library_dirs=['/home/kesler/lib'],
                    libraries =  ['blitz', 'xml2', 'hdf5'])

setup (name = 'IO',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1],
       py_modules=['IO'])
