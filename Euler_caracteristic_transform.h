#include <iostream>
#include <stdexcept>

#ifndef EULER_CARACTERISTIC_TRANSFORM_H_
#define EULER_CARACTERISTIC_TRANSFORM_H_

class Euler_caracteristic_transform{
public:
	std::vector<double> T;
	std::vector<double> transform_values;

	// Passer par reference
	Euler_caracteristic_transform(std::vector<double> sorted_scalar_products, std::vector<double> values){
		T = sorted_scalar_products;
		transform_values = values;
	}

	// pareil
	std::size_t dichotomie(std::vector<double>& table, double t, std::size_t begin, std::size_t end){
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

	double evaluate(double t){
		if(T[T.size()-1] < t){
			return transform_values[transform_values.size()-1];
		}else{
			return transform_values[dichotomie(T, t, 0, T.size())];
		}
	}

	std::vector<double> vectorize(double a, double b, int n){
		if(b <= a){
			std::cout << "[" << a << ";" << b << "] is not a valid interval" << std::endl;
			throw std::invalid_argument("Invalid interval");
		}
		if(n <= 0){
			std::cout << "Cannot vectorize for " << n << " points" << std::endl;
			throw std::invalid_argument("Invalid number of points");
		}
		std::vector<double> res;
		if(n==1){
			res.push_back(evaluate(a));
			return res;
		}

		double step = (b-a) / (n-1);
 		int tid = 0;
		int i = 0;
		while(i<n){
			double t = a + step*i;
			if(tid < (int)T.size()-1){
				if(T[tid+1] > t){
					res.push_back(transform_values[tid]);
					i++;
				}else{
					tid++;
				}
			}else{
				res.push_back(transform_values[tid]);
				i++;
			}
		}
		return res;
	}

	std::vector<std::vector<double>> get_attributes(){
		std::vector<std::vector<double>> v;
		v.push_back(T);
		v.push_back(transform_values);
		return v;
	}
};

#endif //EULER_CARACTERISTIC_TRANSFORM_H_