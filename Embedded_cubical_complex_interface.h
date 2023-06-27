#ifndef INCLUDE_EMBEDDED_CUBICAL_COMPLEX_INTERFACE_H_
#define INCLUDE_EMBEDDED_CUBICAL_COMPLEX_INTERFACE_H_

#include "Embedded_cubical_complex.h"

#include <iostream>
#include <vector>
#include <string>

template<typename Embedded_cubical_complex_options = Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double>>
class Embedded_cubical_complex_interface : public Embedded_cubical_complex<Embedded_cubical_complex_options>{
    public:
    Embedded_cubical_complex_interface(const std::vector<unsigned>& dimensions,
                        const std::vector<double>& cells, bool input_top_cells = true)
    : Embedded_cubical_complex<Embedded_cubical_complex_options>(dimensions, cells, input_top_cells) {
    }
};

#endif