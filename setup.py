from setuptools import setup, Extension
from Cython.Build import cythonize

gudhi_path = "/usr/include"

setup(ext_modules = cythonize(Extension(
        "embedded_cubical_complex",
        sources=["embedded_cubical_complex.pyx"],
        language="c++",
        include_dirs = [gudhi_path]
)))