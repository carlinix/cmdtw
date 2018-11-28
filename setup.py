from distutils.core import setup
from Cython.Build import cythonize
#from setuptools import setup, find_packages, Extension
#from Cython.Distutils import build_ext
from codecs import open
import numpy

setup(
    include_dirs=[numpy.get_include(), 'cdtw/'],
    ext_modules = cythonize("cymdtw.pyx")
)
