#include <time.h>
#include "../../Embedded_cubical_complex.h"

void print_error(){
	std::cout << "This program is made to benchmark the radon transform algorithm.\n"
    <<"This program will create a few complexe and some vectors and compute the radon transform of the complexe in the vectors directions.\n"
    <<"./main n s_0 s_1 ... s_{n-1} r k m \n"
    <<"Arguments : \n"
    <<"- n : the dimension of the complexe.\n"
    <<"- s_0,...,s_{n-1} : the sizes of your complex.\n"
    <<"- r : the complexes will have values in [|0,r_1|].\n"
    <<"- k : the vectors will have their coefficients in [-k,k].\n"
    <<"- m : the number of vectors to create for each complex.\n";
}

std::vector<double> create_data(std::vector<unsigned> sizes, int range){
    std::vector<double> data;
    unsigned num_cells = 1;
    for(size_t i=0; i<sizes.size(); i++){
        num_cells *= sizes[i];
    }
    for(unsigned i=0; i<num_cells; i++){
        data.push_back(std::rand() % range);
    }
    return data;
}

std::vector<double> random_vector(int dimension, double range){
    std::vector<double> vect;
    for(int i=0; i<dimension; i++){
        vect.push_back(((double)std::rand() / RAND_MAX - 0.5) * 2 * range);
    }
    return vect;
}

void recap_parameters(int dimension, std::vector<unsigned> sizes, int range, int vect_count, int vect_range){
	std::cout << "Dimension     		  : " << dimension << "\n";
	std::cout << "Complex sizes 		  : ";
	print_vector(sizes);
	std::cout << "Complex filration range : " << range << "\n";
	std::cout << "Number of vectors       : " << vect_count << "\n";
	std::cout << "Vector range            : " << vect_range << "\n\n";
}

int main(int argc, char** argv){
	if(argc > 2){
		int dimension = std::stoi(argv[1]);
		if(argc > dimension + 2){
			std::vector<unsigned> sizes;
			for(int i=0; i<dimension; i++){
				sizes.push_back(std::stoi(argv[i+2]));
			}
			if(argc > dimension + 4){
				int range = std::stoi(argv[dimension+2]);
				int vect_count = std::stoi(argv[dimension+3]);
				int vect_range = std::stoi(argv[dimension+4]);
				std::vector<std::vector<double>> vectors;
				for(int i=0; i<vect_count; i++){
					vectors.push_back(random_vector(dimension,vect_range));
				}

				//recap_parameters(dimension,sizes,range,vect_count,vect_range);
				Embedded_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>> cplx(sizes, create_data(sizes,range));

				auto t0 = std::chrono::high_resolution_clock::now();
	            cplx.init_radon_transform();
	            auto t1 = std::chrono::high_resolution_clock::now();
	            double tot_precalc = std::chrono::duration<double, std::milli>(t1-t0).count();

				double tot = 0;
				for(int i=0; i<vect_count; i++){
					auto t0 = std::chrono::high_resolution_clock::now();
	            	cplx.compute_radon_transform(vectors[i]);
	            	auto t1 = std::chrono::high_resolution_clock::now();
	            	tot += std::chrono::duration<double, std::milli>(t1-t0).count();
				}

				/*
				std::cout << "Total pre-calculation time : " << tot_precalc << " ms\n";

				std::cout << "Total computation time     : " << tot << " ms\n";
				std::cout << "Average computation time   : " << tot / vect_count << " ms\n";
				*/

				std::cout << tot_precalc << "," << tot/vect_count << "\n";
			}else{
				print_error();
			}
		}else{
			print_error();
		}
	}else{
		print_error();
	}
}