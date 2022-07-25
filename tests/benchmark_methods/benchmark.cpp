#include <time.h>
#include "embedded_complex.hpp"


void message_invalid_arg_number(){
    std::cout << "This program is made to benchmark the different methods to compute hybrid transforms.\n"
    <<"This program will create a fiwed number of complexes and vectors and compute the hybrid transform with kernel t -> exp(t) of the complexes with the generated vectors.\n"
    <<"./main b n s_0 s_1 ... s_{n-1} r p k m \n"
    <<"Arguments : \n"
    <<"- b : if zero, runs all the methods on one core, if 1 runs the critical points method with 1,2,3... cores up to the number of cores that you have.\n"
    <<"- n : the dimension of the complexes.\n"
    <<"- s_0,...,s_{n-1} : the sizes of your complex.\n"
    <<"- r : the complexes will have values in [|0,r_1|].\n"
    <<"- p : the number of complexes to create. \n"
    <<"- k : the vectors will have their coefficients in [-k,k].\n"
    <<"- m : the number of vectors to create for each complex.\n"
    <<"Example : ./main 0 2 10 10 5 100 3 10000 will create 100 10*10 complex with filtration values in [|0,4|] and compute the transform for 10 000 vectors with coordinates between -3 and 3 on each complex\n";
}

void recap_params(std::vector<unsigned> sizes, int complex_num, int iterations, double vector_range, int filtration_range){
    std::cout << "----Benchmark parameters----\n";
    std::cout << "Number of complexes        : " << complex_num << "\n";
    std::cout << "Complexes sizes            : ";
    print_vector(sizes);
    std::cout << "Complexes filtration range : " << filtration_range << "\n";
    std::cout << "Number of vectors          : " << iterations << "\n";
    std::cout << "Vectors range              : " << vector_range << "\n\n";
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

void benchmark_multithread(std::vector<unsigned> sizes, std::vector<std::vector<double>> data, int iterations, double vector_range, int filtration_range){
    recap_params(sizes, data.size(), iterations, vector_range, filtration_range);
    std::vector<double> results_calc;
    std::vector<double> results_pre_calc;

    int dimension = sizes.size();
    std::vector<std::vector<double>> vect_list;
    for(int i=0; i<iterations; i++){
        vect_list.push_back(random_vector(dimension, vector_range));
    }

    for(int i=1; i<=std::thread::hardware_concurrency(); i++){
        double tot0 = 0;
        double tot1 = 0;
        for(int j=0; j<data.size(); j++){
            Embedded_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>> cplx(sizes, data[j]);
            auto t0 = std::chrono::system_clock::now();
            cplx.compute_non_singular_critical_vertices(i);
            auto t1 = std::chrono::system_clock::now();
            cplx.compute_hybrid_transform(&std::exp, vect_list, i);
            auto t2 = std::chrono::high_resolution_clock::now();
            tot0 += std::chrono::duration<double, std::milli>(t1-t0).count();
            tot1 += std::chrono::duration<double, std::milli>(t2-t1).count();
        }
        results_pre_calc.push_back(tot0/data.size());
        results_calc.push_back(tot1/(iterations*data.size()));

        std::cout << "Number of cores used                 : " << i << "\n";
        std::cout << "Time took by pre-calculation         : " << tot0 << "ms\n";
        std::cout << "Average time took by pre-calculation : " << tot0/data.size() << "ms\n";
        std::cout << "Time took by calculation             : " << tot1 << "ms\n";
        std::cout << "Average time took by calculation     : " << tot1/iterations << "ms\n";
        std::cout << "-----------------------------------------\n\n"; 

    }
    std::cout << "Recap : \n";
    std::cout << "Pre-calculations durations : ";
    print_vector(results_pre_calc);
    std::cout << "Average calculations durations : ";
    print_vector(results_calc);
}

void benchmark(std::vector<unsigned> sizes, std::vector<std::vector<double>> data, int iterations, double vector_range, int filtration_range){
    recap_params(sizes, data.size(), iterations, vector_range, filtration_range);
    int dimension = sizes.size();
    std::vector<std::vector<double>> vect_list;
    for(int i=0; i<iterations; i++){
        vect_list.push_back(random_vector(dimension, vector_range));
    }
    
    double tot0 = 0;
    double tot1 = 0;
    double tot2 = 0;
    double tot3 = 0;
    double tot4 = 0;
    double tot5 = 0;

    double diff_simple_arr = 0;
    double diff_simple_crit = 0;

    for(int i=0; i<data.size(); i++){
        Embedded_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>> cplx(sizes, data[i]);

        std::vector<double> simple_res;
        std::vector<double> arrangement_res;
        std::vector<double> crit_res_mono;

        auto t0 = std::chrono::high_resolution_clock::now();
        for(int i=0; i<vect_list.size(); i++){
            simple_res.push_back(cplx.compute_hybrid_transform_simple(&std::exp, vect_list[i]));
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        cplx.compute_arrangement();
        auto t2 = std::chrono::system_clock::now();
        for(int i=0; i<vect_list.size(); i++){
            arrangement_res.push_back(cplx.compute_hybrid_transform_arrangement(&std::exp, vect_list[i]));
        }
        auto t3 = std::chrono::system_clock::now();
        cplx.compute_non_singular_critical_vertices(1);
        auto t4 = std::chrono::system_clock::now();
        crit_res_mono = cplx.compute_hybrid_transform(&std::exp, vect_list, 1);
        auto t5 = std::chrono::high_resolution_clock::now();

        for(int i=0; i<simple_res.size(); i++){ //IL FAUT CALCULER UN TRUC DU TYPE ECART TYPE SINON ON VA PAS Y ARRIVER....
            diff_simple_arr += std::abs(simple_res[i] - arrangement_res[i]);
            diff_simple_crit += std::abs(simple_res[i] - crit_res_mono[i]);
        }

        tot0 += std::chrono::duration<double, std::milli>(t1-t0).count();
        tot1 += std::chrono::duration<double, std::milli>(t2-t1).count();
        tot2 += std::chrono::duration<double, std::milli>(t3-t2).count();
        tot3 += std::chrono::duration<double, std::milli>(t4-t3).count();
        tot4 += std::chrono::duration<double, std::milli>(t5-t4).count();
    }

    std::cout << "-----------------------------------------\n";
    std::cout << "Total computation times : \n";
    std::cout << "Simple method                   : " << tot0 << "ms\n\n";

    std::cout << "Arrangements pre-calculation    : " << tot1 << "ms\n";
    std::cout << "Arrangements calculation        : " << tot2 << "ms\n";
    std::cout << "Arrangements method             : " << tot2+tot1 << "ms\n\n";


    std::cout << "Critical points pre-calculation : " << tot3 << "ms\n";
    std::cout << "Critical points calculation     : " << tot4 << "ms\n";
    std::cout << "Critical points method          : " << tot4+tot3 << "ms\n";

    std::cout << "-----------------------------------------\n\n"; 

    std::cout << "Average computation times : \n";
    std::cout << "Simple method                   : " << (tot0)/iterations << "ms\n";
    std::cout << "Arrangements method             : " << (tot2+tot1)/iterations << "ms\n";
    std::cout << "Critical points method          : " << (tot4+tot3)/iterations << "ms\n";

    std::cout << "-----------------------------------------\n\n";
    std::cout << "Average difference simple/arrangements : " << diff_simple_arr/iterations << "\n";
    std::cout << "Average difference simple/critical-points : " << diff_simple_crit/iterations << "\n";
    std::cout << "-----------------------------------------\n\n";
}   


int main(int argc, char** argv){
    if(argc > 2 && (std::stoi(argv[1]) == 0 || std::stoi(argv[1]) == 1)){
        int dimension = std::stoi(argv[2]);
        if(argc == dimension + 7){
            std::vector<unsigned> sizes;
            for(int i=0; i<dimension; i++){
                sizes.push_back(std::stoi(argv[i+3]));
            }
            int range = std::stoi(argv[dimension+3]);

            std::vector<std::vector<double>> data;
            for(int i=0; i<std::stoi(argv[dimension+4]); i++){
                data.push_back(create_data(sizes, range));
            }
        
            double vector_range = std::stod(argv[dimension+5]);
            int iterations = std::stoi(argv[dimension+6]);
            if(std::stoi(argv[1]) == 0){
                benchmark(sizes, data, iterations, vector_range, range);
            }else{
                benchmark_multithread(sizes, data, iterations, vector_range, range);
            }
        }else{
            message_invalid_arg_number();
        }
    }else{
        message_invalid_arg_number();          
    }
    return 0;
}

