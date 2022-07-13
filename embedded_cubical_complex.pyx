# distutils: language = c++

from libcpp.vector cimport vector
from libc.math cimport exp, sin, cos

import numpy as np

cdef extern from "Radon_transform.h":
    cdef cppclass Radon_transform_base "Radon_transform":
        Radon_transform_base(vector[double] T, vector[double] Values, vector[double] singular_T, vector[double] singular_values) nogil
        double evaluate(double t) nogil
        vector[vector[double]] get_attributes() nogil

cdef extern from "Embedded_cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Embedded_cubical_complex_base_interface "Embedded_cubical_complex_interface<>":
        Embedded_cubical_complex_base_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells) nogil
        void compute_critical_vertices(int num_jobs) nogil
        void print_filtration() nogil
        vector[double] compute_hybrid_transform(int kernel_num, vector[vector[double]] directions_list, int num_jobs) nogil
        vector[vector[double]] compute_radon_transform_python(vector[double] direction) nogil


cdef class RadonTransform:
    cdef Radon_transform_base * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, T, Values, singular_T, singular_values):
        """RadonTransform constructor.
        DO NOT USE THIS TO CONSTRUCT RADON TRANSFORM, USE EmbeddedComplex.compute_radon_transform(direction) INSTEAD 
        """
    
    # The real cython constructor
    def __cinit__(self, T, Values, singular_T, singular_values):
        self.this_ptr = new Radon_transform_base(T, Values, singular_T, singular_values)

    def evaluate(self, t):
        return self.this_ptr.evaluate(t)

    def get_attributes(self):
        return self.this_ptr.get_attributes()

cdef class EmbeddedComplex:

    cdef Embedded_cubical_complex_base_interface * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, top_dimensional_cells=None, dimensions=None, int num_jobs = 0):
        """EmbeddedComplex constructor.
        :param top_dimensional_cells: The filtration values of the top dimensional cells of the cubical complex.
        :type top_dimensional_cells: list of double or numpy ndarray of double 
        :param dimensions: The sizes of the embedded complex. Needed only if top_dimensional_cells is a list of double.
        :type dimensions: list of int
        :param num_jobs:The number of threads to use to compute the critical points of the cubical complex.
        :type num_jobs: integer
        """
    
    # The real cython constructor
    def __cinit__(self, top_dimensional_cells=None, dimensions=None, int num_jobs=0):
        if(top_dimensional_cells is not None):
            if(type(top_dimensional_cells) is np.ndarray):
                self._construct_from_cells(top_dimensional_cells.shape, top_dimensional_cells.ravel(), num_jobs)
            elif(dimensions is not None):
                self._construct_from_cells(dimensions, top_dimensional_cells, num_jobs)            

    def __dealloc__(self):
        if self.this_ptr != NULL:
            del self.this_ptr

    def _find_kernel(self, kernel_name):
        if(kernel_name == "exp"):
            return 0
        elif(kernel_name == "cos"):
            return 1
        elif(kernel_name == "sin"):
            return 2
        return -1       

    def _construct_from_cells(self, vector[unsigned] dimensions, vector[double] top_dimensional_cells, int num_jobs):
        with nogil:
            self.this_ptr = new Embedded_cubical_complex_base_interface(dimensions, top_dimensional_cells)
            self.this_ptr.compute_critical_vertices(num_jobs)

    def compute_hybrid_transform(self, kernel_name, vector[vector[double]] directions, int num_jobs = -1):
        kernel_num = self._find_kernel(kernel_name)
        return self.this_ptr.compute_hybrid_transform(kernel_num, directions, num_jobs)

    def compute_radon_transform(self, vector[double] direction):
        tmp = self.this_ptr.compute_radon_transform_python(direction)
        if len(tmp) == 4:
            radon = RadonTransform(tmp[0],tmp[1],tmp[2],tmp[3])
            return radon

    def print_filtration(self):
        self.this_ptr.print_filtration()