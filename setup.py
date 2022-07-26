from setuptools import setup, Extension
from Cython.Build import cythonize
from numpy import get_include as numpy_get_include

# You have to put the path to the required libraries as a list of strings
libs = ["/usr/include"]

setup(ext_modules = cythonize(Extension(
        "embedded_cubical_complex",
        sources=["embedded_cubical_complex.pyx"],
        language="c++",
        include_dirs = libs,
        extra_compile_args=["-std=c++2a"]
)))