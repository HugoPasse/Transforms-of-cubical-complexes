//#include "Embedded_cubical_complex.h"
#include <vector>

#include <stdexcept>
#include <limits>


#ifndef RADON_TRANSFORM_H_
#define RADON_TRANSFORM_H_

class Radon_transform{
    public :
        std::vector<double> T;
        std::vector<double> Values;

        std::vector<double> singular_T;
        std::vector<double> singular_values;

        std::size_t dichotomie(std::vector<double> table, double t, std::size_t begin, std::size_t end){
            if(end - begin <= 1){
                return begin;
            }else{
                std::size_t index = (begin + end) / 2;
                if(table[index] > t){
                    return dichotomie(table, t, begin, index);
                }else{
                    return dichotomie(table, t, index, end);
                }
            }
        }

        double get_singular_critical_value(double t, bool* is_value_critical){
            std::size_t index = dichotomie(singular_T, t, 0, singular_T.size());
            if(std::abs(singular_T[index] - t) <= std::numeric_limits<double>::epsilon()){
                *is_value_critical = true;
                return singular_values[index];
            }else if(index < singular_T.size() - 1 && std::abs(singular_T[index + 1] - t) <= std::numeric_limits<double>::epsilon()){
                *is_value_critical = true;
                return singular_values[index + 1];
            }
            *is_value_critical = false;
            return 0;
        }

        Radon_transform(std::vector<double> _T, std::vector<double> _Values, std::vector<double> _singular_T, std::vector<double> _singular_Values){
            T = _T;
            Values = _Values;
            singular_T = _singular_T;
            singular_values = _singular_Values;
        }
   
        /*
        Radon_transform(Embedded_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>> cplx, std::vector<double> direction){
            int index = cplx.get_vector_index(direction);
            int reverse_multiplicity = -1;

            if(index >= (int)cplx.critical_vertices.size()){
                reverse_multiplicity = -1;
                index = (1 << direction.size()) - 1 - index;
            }

            std::vector<double> scalar_pdt;
            std::vector<double> singular_scalar_pdt;

            std::vector<std::size_t> indices;
            std::vector<std::size_t> singular_indices;
            
            //Computing scalar products with non singular critical vertices
            for(std::size_t i = 0; i < cplx.critical_vertices[index].size(); i++){
                scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), cplx.embedding[cplx.embedding_index[cplx.critical_vertices[index][i]]].begin(), 0.0));
                indices.push_back(i);
            }
            //Computing scalar products with singular critical vertices
            for(std::size_t i = 0; i < cplx.zero_measure_critical_vertices[index].size(); i++){
                singular_scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), cplx.embedding[cplx.embedding_index[cplx.zero_measure_critical_vertices[index][i]]].begin(), 0.0));
                singular_indices.push_back(i);
            }

            std::sort(indices.begin(), indices.end(), [&scalar_pdt](int i, int j) {return scalar_pdt[i] < scalar_pdt[j];});
            std::sort(singular_indices.begin(), singular_indices.end(), [&singular_scalar_pdt](int i, int j) {return singular_scalar_pdt[i] < singular_scalar_pdt[j];});
            
            //Filling T with changing points and Values[i] = Radon(t) for all t \in [T[i],T[i+1]]
            //Last element of values should always be 0
            std::size_t len = 0;
            int euler_car = reverse_multiplicity * cplx.critical_multiplicity[index][indices[0]];
            T.push_back(scalar_pdt[indices[0]]);
            Values.push_back(euler_car);
            for(std::size_t i = 1; i < indices.size(); i++){        
                int crit_mul = reverse_multiplicity * cplx.critical_multiplicity[index][indices[i]];
                euler_car += crit_mul;
                
                if(std::abs(scalar_pdt[indices[i-1]] - scalar_pdt[indices[i]]) <= std::numeric_limits<double>::epsilon()){
                    Values[len] = Values[len] + crit_mul;
                }else{
                    T.push_back(scalar_pdt[indices[i]]);
                    Values.push_back(euler_car);
                    len++;
                }
            }

            if(singular_indices.size() > 0){
                len = 0;
                singular_T.push_back(singular_scalar_pdt[singular_indices[0]]);
                euler_car = Values[dichotomie(T,singular_T[0],0,T.size())];
                singular_values.push_back(euler_car - cplx.zero_measure_critical_multiplicity[index][singular_indices[0]]);

                for(std::size_t i = 1; i < singular_indices.size(); i++){        
                    int crit_mul = cplx.zero_measure_critical_multiplicity[index][singular_indices[i]];

                    if(std::abs(singular_scalar_pdt[indices[i-1]] - singular_scalar_pdt[singular_indices[i]]) <= std::numeric_limits<double>::epsilon()){
                        singular_values[len] = singular_values[len] + crit_mul;
                    }else{
                        singular_T.push_back(singular_scalar_pdt[singular_indices[i]]);
                        len++;
                        singular_values.push_back(Values[dichotomie(T,singular_T[len],0,T.size())] - cplx.zero_measure_critical_multiplicity[index][singular_indices[len]]);
                    }
                }
            }
        }*/

        double evaluate(double t){
            bool is_value_critical = false;
            double val = get_singular_critical_value(t, &is_value_critical);

            if(is_value_critical){
                return val;
            }

            if(t < T[0]){
                return 0;
            }else{
                int index = dichotomie(T,t,0,T.size());
                return Values[index];
            }
        }

        std::vector<std::vector<double>> get_attributes(){
            std::vector<std::vector<double>> result;

            result.push_back(T);
            result.push_back(Values);
            result.push_back(singular_T);
            result.push_back(singular_values);

            return result;
        }

        /*
        void print_exact_radon(){
            std::cout << "T :\n";
                print_vector(T);
                std::cout << "Values :\n";
                print_vector(Values);
                std::cout << "\n";

                std::cout << "singular_T :\n";
                print_vector(singular_T);
                std::cout << "singular_Values :\n";
                print_vector(singular_values);
                std::cout << "\n";
        }*/
};

#endif //define RADON_TRANSFORM_H_