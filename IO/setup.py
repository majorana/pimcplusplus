from distutils.core import setup, Extension

module1 = Extension('IOSection',
                    sources = ['InputOutputPythonWrapper.cc','InputOutput.cc',\
                               'InputOutputHDF5.cc', 'InputOutputASCII.cc',\
                               'InputOutputXML.cc' ],\
                    include_dirs=['/home/common/lib/blitz-0.6-GCC/',\
                                  '/home/common/lib/hdf5/include',
                                  '/usr/include/libxml2'],
                    library_dirs=['/home/common/lib/blitz-0.6-GCC/lib',\
                                  '/home/common/lib/hdf5/lib/'],
                    libraries =  ['blitz', 'xml2', 'hdf5'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])
