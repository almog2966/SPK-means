from setuptools import setup, Extension
module1 = Extension('myspkmeans', sources=['spkmeansmodule.c', 'spkmeans.c'])
module = Extension('kmeans_c', sources=['kmeansmodule.c','kmeans.c'])
setup(name='myspkmeans',
      version='1.0',
      description='spkmeans algorithm',
      ext_modules=[module1,module])