# distutils: language = c++

from libcpp.vector cimport vector
from libc.math cimport exp, sin, cos
from libcpp cimport bool

import numpy as np

cdef extern from "Radon_transform.h":
    cdef cppclass Radon_transform_base "Radon_transform":
        Radon_transform_base(vector[double] T, vector[double] Values, vector[double] singular_T, vector[double] singular_values) nogil
        double evaluate(double t) nogil
        vector[vector[double]] get_attributes() nogil

cdef extern from "Euler_caracteristic_transform.h":
    cdef cppclass ECT_base "Euler_caracteristic_transform":
        ECT_base(vector[double] T, vector[double] sorted_transform_values) nogil
        double evaluate(double t) nogil
        vector[vector[double]] get_attributes() nogil

cdef extern from "Embedded_cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Embedded_cubical_complex_base_interface "Embedded_cubical_complex_interface<>":
        Embedded_cubical_complex_base_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells, bool input_top_cells) nogil
        void impose_lower_star_filtration() nogil
        void impose_upper_star_filtration() nogil
        void impose_lower_star_filtration_from_vertices() nogil

        void preproc_hybrid_transform(int num_jobs) nogil
        void preproc_radon_transform(int num_jobs) nogil
        void preproc_ect(int num_jobs) nogil

        vector[double] compute_hybrid_transform(int kernel_num, vector[vector[double]] directions_list, int num_jobs) nogil
        vector[vector[double]] compute_radon_transform_python(vector[double] direction) nogil
        vector[vector[double]] compute_ect_python(vector[double] direction) nogil
        double compute_euler_caracteristic_of_complex() nogil
        
        int get_vector_index(vector[double] direction) nogil
        
        vector[int] get_ordinary_critical_vertices(int index) nogil
        vector[int] get_ordinary_critical_values(int index) nogil
        vector[int] get_classical_critical_vertices(int index) nogil
        vector[int] get_classical_critical_values(int index) nogil

        vector[double] get_vertex_coordinates(int index) nogil
        
        void print_filtration() nogil
        void print_embedding() nogil

cdef class RadonTransform:
    cdef Radon_transform_base * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, T, Values, singular_T, singular_values):
        """RadonTransform constructor.
        DO NOT USE THIS TO CONSTRUCT THE RADON TRANSFORM, USE EmbeddedComplex.compute_radon_transform(direction) INSTEAD 
        """
    
    # The real cython constructor
    def __cinit__(self, T, Values, singular_T, singular_values):
        self.this_ptr = new Radon_transform_base(T, Values, singular_T, singular_values)

    def evaluate(self, t):
        return self.this_ptr.evaluate(t)

    def get_attributes(self):
        return self.this_ptr.get_attributes()

cdef class EulerCaracteristicTransform:
    cdef ECT_base * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, T, transform_values):
        """RadonTransform constructor.
        DO NOT USE THIS TO CONSTRUCT THE EULER CARACTERISTIC TRANSFORM, USE EmbeddedComplex.compute_euler_caracteristic_transform(direction) INSTEAD 
        """
    
    # The real cython constructor
    def __cinit__(self, T, transform_values):
        self.this_ptr = new ECT_base(T, transform_values)

    def evaluate(self, t):
        return self.this_ptr.evaluate(t)

    def get_attributes(self):
        return self.this_ptr.get_attributes()

cdef class EmbeddedComplex:

    cdef Embedded_cubical_complex_base_interface * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, top_dimensional_cells=None, dimensions=None, input_top_cells=True, int num_jobs = 0):
        """EmbeddedComplex constructor.
        :param top_dimensional_cells: The filtration values of the top dimensional cells of the cubical complex.
        :type top_dimensional_cells: list of double or numpy ndarray of double 
        :param dimensions: The sizes of the embedded complex. Needed only if top_dimensional_cells is a list of double.
        :type dimensions: list of int
        :param num_jobs:The number of threads to use to compute the critical points of the cubical complex.
        :type num_jobs: integer
        """
    
    # The real cython constructor
    def __cinit__(self, top_dimensional_cells=None, dimensions=None, input_top_cells=True, int num_jobs=0):
        if(top_dimensional_cells is not None):
            if(type(top_dimensional_cells) is np.ndarray):
                self._construct_from_cells(top_dimensional_cells.shape[::-1], top_dimensional_cells.ravel(), input_top_cells, num_jobs)
            elif(dimensions is not None):
                self._construct_from_cells(dimensions, top_dimensional_cells, input_top_cells, num_jobs)            

    def impose_lower_star_filtration(self):
        self.this_ptr.impose_lower_star_filtration()

    def impose_upper_star_filtration(self):
        self.this_ptr.impose_upper_star_filtration()

    def impose_lower_star_filtration_from_vertices(self):
        self.this_ptr.impose_lower_star_filtration_from_vertices()

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

    def _construct_from_cells(self, vector[unsigned] dimensions, vector[double] top_dimensional_cells, bool input_top_cells, int num_jobs):
        with nogil:
            self.this_ptr = new Embedded_cubical_complex_base_interface(dimensions, top_dimensional_cells, input_top_cells)

    def preproc_hybrid_transform(self,int num_jobs=0):
        self.this_ptr.preproc_hybrid_transform(num_jobs)
    
    def preproc_radon_transform(self,int num_jobs=0):
        self.this_ptr.preproc_radon_transform(num_jobs)
    
    def preproc_ect(self,int num_jobs=0):
        self.this_ptr.preproc_ect(num_jobs)

    def compute_hybrid_transform(self, kernel, vector[vector[double]] directions, int num_jobs = -1):
        if isinstance(kernel, str):
            kernel_num = self._find_kernel(kernel)
            return np.array(self.this_ptr.compute_hybrid_transform(kernel_num, directions, num_jobs))
        else:
            R = np.zeros(directions.size())
            for k in range(len(directions)):
                d = directions[k]
                S = 0
                index = self.get_vector_index(d)
                crit = self.get_critical_vertices(index)
                mult = self.get_critical_multiplicity(index)
                for i in range(len(crit)):
                    v = crit[i]
                    m = mult[i]
                    S += m * kernel(np.dot(self.get_coordinates(v),d))
                R[k] = S
            return R


    def compute_radon_transform(self, vector[double] direction):
        tmp = self.this_ptr.compute_radon_transform_python(direction)
        radon = RadonTransform(tmp[0],tmp[1],tmp[2],tmp[3])
        return radon

    def compute_euler_caracteristic_transform(self, vector[double] direction):
        tmp = self.this_ptr.compute_ect_python(direction)
        ect = EulerCaracteristicTransform(tmp[0],tmp[1])
        return ect

    def compute_euler_caracteristic_of_complex(self):
        return self.this_ptr.compute_euler_caracteristic_of_complex()

    def get_ordinary_critical_vertices(self, index):
        return self.this_ptr.get_ordinary_critical_vertices(index)

    def get_ordinary_critical_values(self, index):
        return self.this_ptr.get_ordinary_critical_values(index)        
    
    def get_classical_critical_vertices(self, index):
        return self.this_ptr.get_classical_critical_vertices(index)        
    
    def get_classical_critical_values(self, index):
        return self.this_ptr.get_classical_critical_values(index)
    
    def get_coordinates(self, int vertex):
        return self.this_ptr.get_vertex_coordinates(vertex)

    def get_vector_index(self, vector[double] v):
        index = 0
        i = 1
        for c in v:
            if c > 0:
                index += i
            i *= 2
        return index

    def print_filtration(self):
        self.this_ptr.print_filtration()

    def print_embedding(self):
        self.this_ptr.print_embedding()
