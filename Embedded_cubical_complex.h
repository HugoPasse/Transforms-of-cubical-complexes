#include <gudhi/Bitmap_cubical_complex.h>

#include <future>
#include <thread>
#include <chrono>
#include <time.h>

#include "Radon_transform.h"
#include "Euler_caracteristic_transform.h"

#ifndef EMBEDDED_CUBICAL_COMPLEX_H_
#define EMBEDDED_CUBICAL_COMPLEX_H_

template <typename T>
void print_vector(std::vector<T> vect){
    if(vect.size() == 0){
        std::cout << "[]\n";
    }else{
        std::cout << "[" << vect[0];
        for(std::size_t i = 1; i < vect.size(); i++){
            std::cout << ", " << vect[i];
        }
        std::cout << "]\n";
    }
}

template <typename T>
class Embedded_cubical_complex : public Gudhi::cubical_complex::Bitmap_cubical_complex<T>
{
    public:
        typedef std::size_t Simplex_key;
        typedef typename T::filtration_type Filtration_value;
        typedef Simplex_key Simplex_handle;
    
        std::vector<std::vector<double>> embedding;         //Array of the cubical complexe's points' coordinates
        std::vector<int> embedding_index;                   //Array to link vertices index in the cubical complex to their index in the embedding

        std::vector<int> sizes_pdt;                         //Products of the sizes from s_0 up to s_0*s_1*...*s_(n-2)

        // Classical critical points
        int are_cla_crit_pts_computed = 0;
        std::vector<std::vector<int>> cla_crit_pts;
        std::vector<std::vector<int>> cla_crit_val;

        // Ordinary critical points
        int are_ord_crit_pts_computed = 0;
        std::vector<std::vector<int>> ord_crit_pts;
        std::vector<std::vector<int>> ord_crit_val;
        //If we have a vector e = (e_1,...,e_n), its index is the sum of the 2^i were i are the indexes such that e_i >= 0
        //If index_v = index_w they have the same critical points. The critical points are stored in critical_vertices[index]

        //*********************************************//
        // Constructor
        //*********************************************//
        Embedded_cubical_complex(const std::vector<unsigned>& dimensions,
            const std::vector<Filtration_value>& cells, bool input_top_cells = true):Gudhi::cubical_complex::Bitmap_cubical_complex<T>(dimensions, cells, input_top_cells)
            {
                sizes_pdt.push_back(2*this->sizes[0]+1);

                //In this loop we compute the product of the number of cubes in each direction, it optimizes the complexity of the get_coordinates_in_complex function
                for(Simplex_handle i=1; i<this->dimension(); i++){
                    sizes_pdt.push_back(sizes_pdt[i-1]*(2*this->sizes[i]+1));
                }

                initalize_embedding();
                initalize_embedding_index();
                //impose_upper_star_filtration();

            }     
        
        //*********************************************//
        // Specific filtration
        //*********************************************//
        void impose_upper_star_filtration(){
            // std::cout << "Imposing upper star\n\n";
            for(Gudhi::cubical_complex::Bitmap_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>>::Top_dimensional_cells_iterator it = this->top_dimensional_cells_iterator_begin(); it != this->top_dimensional_cells_iterator_end(); ++it){
                std::vector<std::size_t> boundary = this->get_boundary_of_a_cell(*it);
                for(std::size_t i=0; i<boundary.size(); i++){
                    //this->data[boundary[i]] = std::max(this->filtration(boundary[i]),this->filtration(*it));
                    impose_upper_star_filtration_from_simplex(boundary[i],this->filtration(*it));
                }
            }
        }

        void impose_upper_star_filtration_from_simplex(Simplex_handle sh, double top_cell_filtration){
            this->data[sh] = std::max(this->filtration(sh),top_cell_filtration);
            if(this->dimension(sh) > 0){
                std::vector<std::size_t> boundary = this->get_boundary_of_a_cell(sh);
                for(std::size_t i=0; i<boundary.size(); i++){
                    //this->data[boundary[i]] = std::max(this->filtration(boundary[i]),this->filtration(sh));
                    impose_upper_star_filtration_from_simplex(boundary[i],top_cell_filtration);
                }
            }
        }

        //*********************************************//
        // Functions for pretratement
        //*********************************************//

        void initalize_embedding(){
            int m = this->sizes[std::distance(this->sizes.begin(), std::max_element(this->sizes.begin(), this->sizes.end()))];
            for(Simplex_handle i = 0; i < this->num_simplices();i++){
                if(this->dimension(i) == 0){
                    std::vector<int> coords = get_coordinates_in_complex(i);
                    std::vector<double> embedded_coords(coords.size());
                    for(Simplex_handle j = 0; j < coords.size(); j++){
                        embedded_coords[j] = (coords[j] - (int)this->sizes[j]) / 2. / m;
                    }
                    embedding.push_back(embedded_coords);
                }
            }
        }

        void initalize_embedding_index(){
            int index = 0;
            for(Simplex_handle handle = 0; handle < this->num_simplices(); handle++){
                if(this->dimension(handle) == 0){
                    embedding_index.push_back(index);
                    index++;
                }else{
                    embedding_index.push_back(-1);
                }
            }
        }

        void preproc_hybrid_transform(int num_jobs=0){
            compute_ordinary_critical_vertices(num_jobs);
        }

        void preproc_radon_transform(int num_jobs=0){
            compute_classical_critical_vertices(num_jobs);
            compute_ordinary_critical_vertices(num_jobs);
        }

        void preproc_ect(int num_jobs=0){
            compute_classical_critical_vertices(num_jobs);
        }

        //*********************************************//
        //Functions for critical points
        //*********************************************//
        void compute_ordinary_critical_vertices(int num_jobs = 0){
            if(are_ord_crit_pts_computed == 0){
                if(num_jobs > (int)std::thread::hardware_concurrency() || num_jobs <= 0){
                    num_jobs = std::thread::hardware_concurrency();
                }
                num_jobs = std::min(num_jobs, (int)this->sizes[0]+1);

                int dimension = this->dimension();
                int n_simplices = this->num_simplices();

                std::vector<int> direction(dimension,-1);
                int index = 0;
                
                for(int i = 0; i < (1 << dimension); i++){    //Loop on every possible direction
                    std::vector<int> tmp;
                    ord_crit_pts.push_back(tmp);
                    ord_crit_val.push_back(tmp);

                    std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                    std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                    std::vector<std::thread> thread_vector;

                    for(int j = 0; j < num_jobs; j++){
                        std::promise<int> promiseObj;
                        promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                        future_vect.push_back(promise_vect[j].get_future());
                        //Objects to get return values from the thread
                        thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_ordinary_critical_vertices_subroutine, this, std::move(promise_vect[j]), direction, dimension, n_simplices, j, num_jobs));
                        thread_vector[j].detach();  //Thread is now running concurently
                    }
                    int job = 0;
                    std::chrono::seconds zero_sec{0};
                    while(job < num_jobs){
                        //Getting thread response and appending it to the corresponding vectors
                        if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                            std::vector<std::vector<int>> thread_res = future_vect[job].get();
                            ord_crit_pts[index].insert(ord_crit_pts[index].end(), thread_res[0].begin(), thread_res[0].end());
                            ord_crit_val[index].insert(ord_crit_val[index].end(), thread_res[1].begin(), thread_res[1].end());
                            job++;
                        }
                    }

                    // Finding next vector
                    int k=0;
                    while(direction[k] == 1 && k < dimension){
                        direction[k] = -1;
                        index -= (1u << k);
                        k++;
                    }
                    direction[k] = 1;
                    index += (1u << k);
                }
            }
            are_ord_crit_pts_computed = 1;
        }

        void compute_ordinary_critical_vertices_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> coords(dim,0);
            Simplex_key vertex = 2*job_index;

            coords[0] = 2*job_index;

            std::vector<int> subroutine_points;
            std::vector<int> subroutine_variations;

            while(vertex < n_simplices){        //Loop on all vertices
                int euler_car = compute_euler_car_in_direction(vertex, direction, 1) - compute_euler_car_in_direction(vertex, direction, -1);
                if(euler_car != 0){
                    subroutine_points.push_back(vertex);
                    subroutine_variations.push_back(euler_car);
                }
                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if(coords[j] > (int)(2*this->sizes[j]+1)){
                        if(j == 0){
                            coords[0] = 2*job_index;
                        }else{
                            coords[j] = 0;
                        }
                        coords[j+1] = coords[j+1] + 2;
                    }else{
                        break;
                    }
                }
                vertex = get_key_from_coordinates(coords);
            }
            std::vector<std::vector<int>> res;
            res.push_back(subroutine_points);
            res.push_back(subroutine_variations);
            promiseObj.set_value(res);
        }

        void compute_classical_critical_vertices(int num_jobs = 0){
            if(are_cla_crit_pts_computed == 0){
                if(num_jobs > (int)std::thread::hardware_concurrency() || num_jobs <= 0){
                    num_jobs = std::thread::hardware_concurrency();
                }
                num_jobs = std::min(num_jobs, (int)this->sizes[0]+1);

                int dimension = this->dimension();
                int n_simplices = this->num_simplices();

                std::vector<int> direction(dimension,-1);
                int index = 0;
                
                for(int i = 0; i < (1 << dimension); i++){    //Loop on every possible direction
                    std::vector<int> tmp;
                    cla_crit_pts.push_back(tmp);
                    cla_crit_val.push_back(tmp);

                    std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                    std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                    std::vector<std::thread> thread_vector;

                    for(int j = 0; j < num_jobs; j++){
                        std::promise<int> promiseObj;
                        promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                        future_vect.push_back(promise_vect[j].get_future());
                        //Objects to get return values from the thread
                        thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_classical_critical_vertices_subroutine, this, std::move(promise_vect[j]), direction, dimension, n_simplices, j, num_jobs));
                        thread_vector[j].detach();  //Thread is now running concurently
                    }
                    int job = 0;
                    std::chrono::seconds zero_sec{0};
                    while(job < num_jobs){
                        //Getting thread response and appending it to the corresponding vectors
                        if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                            std::vector<std::vector<int>> thread_res = future_vect[job].get();
                            cla_crit_pts[index].insert(cla_crit_pts[index].end(), thread_res[0].begin(), thread_res[0].end());
                            cla_crit_val[index].insert(cla_crit_val[index].end(), thread_res[1].begin(), thread_res[1].end());
                            job++;
                        }
                    }

                    // Finding next vector
                    int k=0;
                    while(direction[k] == 1 && k < dimension){
                        direction[k] = -1;
                        index -= (1u << k);
                        k++;
                    }
                    direction[k] = 1;
                    index += (1u << k);
                }
            }
            are_cla_crit_pts_computed = 1;
        }

        void compute_classical_critical_vertices_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> coords(dim,0);
            Simplex_key vertex = 2*job_index;

            coords[0] = 2*job_index;

            std::vector<int> subroutine_points;
            std::vector<int> subroutine_variations;

            while(vertex < n_simplices){        //Loop on all vertices
                int euler_car = this->filtration(vertex) + compute_euler_car_in_direction(vertex, direction, -1);
                if(euler_car != 0){
                    subroutine_points.push_back(vertex);
                    subroutine_variations.push_back(euler_car);
                }
                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if(coords[j] > (int)(2*this->sizes[j]+1)){
                        if(j == 0){
                            coords[0] = 2*job_index;
                        }else{
                            coords[j] = 0;
                        }
                        coords[j+1] = coords[j+1] + 2;
                    }else{
                        break;
                    }
                }
                vertex = get_key_from_coordinates(coords);
            }
            std::vector<std::vector<int>> res;
            res.push_back(subroutine_points);
            res.push_back(subroutine_variations);
            promiseObj.set_value(res);
        }

        //*********************************************//
        //Printing functions
        //*********************************************//

        void print_embedding(){
            std::cout << "[";
            for(Simplex_handle i = 0; i < embedding_index.size(); i++){
                if(embedding_index[i] != -1){
                    print_vector(embedding[embedding_index[i]]);
                }
            }
            std::cout << "]\n";
        }

        void print_filtration(){
            std::clog << "Filtration : \n[";
            for(int i=2*this->sizes[1]; i>=0; i--){
                for(int j=0; j<(int)(2*this->sizes[0]+1); j++){
                    std::clog << this->filtration(j+i*(2*this->sizes[0]+1)) << ", ";
                }
                std::clog << "]\n";
            }
            std::clog << "]\n";
        }

        //*********************************************//
        //Those are arithmetic functions to find the index of the vertex within the vertices
        //*********************************************//

        //This function gives the coordinates of the cell given by key
        //The coordinates are the cartesians' ones where the complex is seen as a cartesian coordinate system
        std::vector<int> get_coordinates_in_complex(Simplex_key key){

            int n = (int)this->dimension();
            std::vector<int> coordinates;      //Coordinates of the face indexed by key to return
            
            for(int i = n-1; i > 0; i--){
                coordinates.insert(coordinates.begin(),key / sizes_pdt[i-1]);
                key = key - coordinates[0]*sizes_pdt[i-1];
            }

            coordinates.insert(coordinates.begin(),key);

            return coordinates;
        }

        //This functions gives you the key of the vertex given the coordinates of it
        //The opposite operation than the previous function
        Simplex_key get_key_from_coordinates(std::vector<int> coordinates){
            Simplex_key key = 0;
            for(int i = this->dimension()-1; i >= 0; i--){
                key = key*(2*this->sizes[i] + 1) + coordinates[i];
            }
            return key;
        }

        //This function returns a vector with the keys of the vertices of the cell given by key
        std::vector<int> get_cell_vertices(Simplex_key key){
            std::vector<int> cell_coordinates = get_coordinates_in_complex(key);    //Get the coordinates of cel indexed by key
            int n = cell_coordinates.size();

            std::vector<int> odd_coordinates;
            std::vector<int> cell_vertex(n);                //The first identified vertex 
            int n_odd_coords = 0;                           //Gives us the dimension of the cell, as well as the number of vertices in it
            
            for(int i = 0; i < n; i++){                         //Computing the number and indexes of odd coordinates in vector cell_coordinates and finding the first vertex
                if(cell_coordinates[i]%2 == 1){
                    odd_coordinates.push_back(i);
                    cell_vertex[i] = cell_coordinates[i]-1;
                    n_odd_coords++;
                }else{
                    cell_vertex[i] = cell_coordinates[i];
                }
            }

            std::vector<int> cell_vertices;

            if(n_odd_coords == 0){  //If key is a vertex we return key
                cell_vertices.push_back(key);
                return cell_vertices;  
            }

            std::vector<int> tmp_vect(n_odd_coords);

            cell_vertices.push_back(get_key_from_coordinates(cell_vertex));

            for(int i = 0; i < (1 << n_odd_coords)-1;i++){
                for(int j = 0; j < n_odd_coords;j++){                   //We use a binary counter on the odd cordinates of the simplex to compute all of the vertices coordinates
                    if(tmp_vect[j] == 1){                           //Cell is vertex if and only if it all of its coordinates are even
                        tmp_vect[j] = 0;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] - 2;
                    }else{
                        tmp_vect[j] = 1;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] + 2;
                        break;
                    }
                }
                
                cell_vertices.push_back(get_key_from_coordinates(cell_vertex));     //Adding the key of the vertex we found in cell_vertices
            }

            return cell_vertices;
        }

        //Given a direction vector e, return the index of the subvector that contains critical points in direction e
        //Maybe these functions must be united in a template one
        template <typename VT>
        int get_vector_index(std::vector<VT> e){
            int index = 0;
            int muliplier = 1;

            for(int i = 0; i < (int)e.size(); i++){
                if(e[i] >= 0){
                    index += muliplier;
                }
                muliplier = muliplier*2;
            }
            return index;
        }

        //Check if given coordinates are in complex, useful when computing the adjacent cells indexes to compute the euler caracteristic
        int are_coordinates_in_complex(std::vector<int> coordinates){
            std::size_t coord_dim = coordinates.size();
            if(coord_dim != this->sizes.size()){
                return 0;
            }

            for(std::size_t i = 0; i < coord_dim; i++){
                if(coordinates[i] < 0 || coordinates[i] > 2*((int)this->sizes[i])){
                    return 0;
                }
            }
            return 1;
        }

        //Return euler car with multiplicity of the intersection between the cells in the neigbourhood of the vertex and the hyperplane orthogonal to direction in the neighbourhood of the vertex
        int compute_euler_car_in_direction(Simplex_handle vertex, std::vector<int> direction, int reverse_vector){
            int euler_car = 0;
            int dim = direction.size();

            std::vector<int> coordinates = get_coordinates_in_complex(vertex); //This vector will successively take the coordinates of adjacent cells involved in the Euler's caracteristic calculation
            
            std::vector<int> tmp(dim);  //This vector will help us to find all the adjacent cells involved in calculations
            int simplex_dim_sign = 1;
            
            for(int i = 0; i < (1 << dim)-1; i++){  //Looping on all of the adjacent cell in direction
                //Finding the next adjacent cell 
                int k = 0;
                while(tmp[k] == 1){             
                    tmp[k] = 0;
                    coordinates[k] -= reverse_vector * direction[k];
                    simplex_dim_sign *= -1;     //Gives us the sign to put in front (e.g : edges must be taken positively were faces must be taken negatively (it is the oppostie to the formulae due to the intersection with the hyperplane that transforms n dimensional cells in n-1 dimensional cells))
                    k++;
                }
                coordinates[k] += reverse_vector * direction[k];

                if(k < dim){
                    tmp[k] = 1;
                    simplex_dim_sign *= -1;
                }

                if(are_coordinates_in_complex(coordinates) == 1){   //If the cell exists, adding a term to the caracteristic
                    Simplex_key key = get_key_from_coordinates(coordinates);
                    euler_car += this->filtration(key) * simplex_dim_sign;
                }
            }
            return euler_car;
        }

        //*********************************************//
        //Functions to compute hybrid transform
        //*********************************************//
        
        //Compute hybrid transform of the complex in direction e, kernel is the antiderivative of the real kernel of the transform
        double compute_hybrid_transform(double (*kernel)(double), std::vector<double> &e){
            if(are_ord_crit_pts_computed == 0){
                compute_ordinary_critical_vertices();
            }

            int index = get_vector_index(e);
            double sum = 0.0;

            for(std::size_t i = 0; i < ord_crit_pts[index].size(); i++){   //Looping on critical vertices
                sum -=  ord_crit_val[index][i] * kernel(std::inner_product(e.begin(),e.end(),embedding[embedding_index[ord_crit_pts[index][i]]].begin(),0.0));
            }
            return sum;
        }

        //An overload of previous function to support multithreading
        std::vector<double> compute_hybrid_transform(double (*kernel)(double), std::vector<std::vector<double>> &vect_list, unsigned int num_threads = -1){
            if(are_ord_crit_pts_computed == 0){
                compute_ordinary_critical_vertices();
            }

            std::vector<double> results;
            std::chrono::seconds zero_sec{0};

            std::size_t num_vectors = vect_list.size();
            
            if(num_threads > std::thread::hardware_concurrency() || num_threads <= 0){
                num_threads = std::thread::hardware_concurrency();
            }

            int step = (int)num_vectors / (int)num_threads;
            
            std::vector<std::promise<std::vector<double>>> promise_vect;
            std::vector<std::future<std::vector<double>>> future_vect;
            std::vector<std::thread> thread_vector;
            
            for(std::size_t i = 0; i < num_threads; i++){   //We create threads
            
                promise_vect.push_back(std::promise<std::vector<double>>());
                future_vect.push_back(promise_vect[i].get_future());
                //Objects to get return values from the thread

                int begin = i*step;
                int end = (i+1)*step;

                if(i == num_threads-1){
                    end = (int)num_vectors;
                }
                
                thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_hybrid_transform_subvector, this, std::move(promise_vect[i]), kernel, std::ref(vect_list), begin, end));
                thread_vector[i].detach();  //Thread is now running concurently
            }

            int b = 1;
            while(b){
                b = 0;
                for(std::size_t i = 0; i < num_threads; i++){       //Waiting for all the threads to finish their job
                    if(future_vect[i].wait_for(zero_sec) != std::future_status::ready){
                        b = 1;
                        break;
                    }
                }
            }
            
            for(std::size_t i = 0; i < num_threads; i++){   //Merging answers in one vector
                std::vector<double> thread_res = future_vect[i].get();
                results.insert(results.end(),thread_res.begin(),thread_res.end());
            }
            
            return results;
        }

        //This overload is made for python wrapping, replacing a pointer to a function by an integer
        std::vector<double> compute_hybrid_transform(int kernel_number, std::vector<std::vector<double>> &vect_list, unsigned int num_threads = -1){
            double (*kernel)(double);
            switch(kernel_number){
                case 0:
                    kernel = &std::exp;
                    break;
                case 1:
                    kernel = &std::cos;
                    break;
                case 2:
                    kernel = &std::sin;
                    break;
                default:
                    throw("Unknown kernel number");
            }
            return compute_hybrid_transform(kernel, vect_list, num_threads);
        }

        //Computing multiple transforms on one kernel, used by previous function, each thread ran by 'compute_hybrid_transform' is going to run an instance of this function
        void compute_hybrid_transform_subvector(std::promise<std::vector<double>> promiseObj, double (*kernel)(double), std::vector<std::vector<double>> &vect_list, std::size_t begin_index, std::size_t end_index){
            std::vector<double> results;
            for(std::size_t i = begin_index; i < end_index; i++){
                results.push_back(compute_hybrid_transform(kernel,vect_list[i]));
            }
            promiseObj.set_value(results);
        }

        //*********************************************//
        //Functions to compute radon transform
        //*********************************************//
        
        Radon_transform compute_radon_transform(std::vector<double> &direction){
            if(are_ord_crit_pts_computed == 0){
                compute_ordinary_critical_vertices();
            }
            if(are_cla_crit_pts_computed == 0){
                compute_classical_critical_vertices();
            }
            std::vector<double> _T;
            std::vector<double> _Values;

            _T.push_back(-std::numeric_limits<double>::infinity());
            _Values.push_back(0);

            std::vector<double> _singular_T;
            std::vector<double> _singular_values;

            int index = get_vector_index(direction);

            std::vector<double> scalar_pdt;
            std::vector<double> singular_scalar_pdt;

            std::vector<std::size_t> indices;
            std::vector<std::size_t> singular_indices;
            
            //Computing scalar products with non singular critical vertices
            for(std::size_t i = 0; i < ord_crit_pts[index].size(); i++){
                scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[ord_crit_pts[index][i]]].begin(), 0.0));
                indices.push_back(i);
            }
            //Computing scalar products with singular critical vertices
            for(std::size_t i = 0; i < cla_crit_pts[index].size(); i++){
                singular_scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[cla_crit_pts[index][i]]].begin(), 0.0));
                singular_indices.push_back(i);
            }

            //We sort the lists of indices because we want to sort critical points by scalar product.
            std::sort(indices.begin(), indices.end(), [&scalar_pdt](int i, int j) {return scalar_pdt[i] < scalar_pdt[j];});
            std::sort(singular_indices.begin(), singular_indices.end(), [&singular_scalar_pdt](int i, int j) {return singular_scalar_pdt[i] < singular_scalar_pdt[j];});
            //Filling T with changing points and Values[i] = Radon(t) for all t \in [T[i],T[i+1]]
            //Last element of values should always be 0
            
            std::size_t len = 0;
            int euler_car = 0;
            //We compute the euler caracteristic values for intervals
            if(indices.size() > 0){ 
                euler_car = ord_crit_val[index][indices[0]];
                _T.push_back(scalar_pdt[indices[0]]);
                _Values.push_back(euler_car);
                for(std::size_t i = 1; i < indices.size(); i++){        
                    int crit_mul = ord_crit_val[index][indices[i]];
                    euler_car += crit_mul;
                    
                    if(std::abs(scalar_pdt[indices[i-1]] - scalar_pdt[indices[i]]) <= std::numeric_limits<double>::epsilon()){
                        _Values[len] = _Values[len] + crit_mul;
                    }else{
                        _T.push_back(scalar_pdt[indices[i]]);
                        _Values.push_back(euler_car);
                        len++;
                    }
                }
            }

            //We compute the euler caracteristic at singular points
            if(singular_indices.size() > 0){
                len = 0;
                _singular_T.push_back(singular_scalar_pdt[singular_indices[0]]);
                _singular_values.push_back(cla_crit_val[index][singular_indices[0]]);
                for(std::size_t i = 1; i < singular_indices.size(); i++){
                    int crit_mul = cla_crit_val[index][singular_indices[i]];
                    if(std::abs(singular_scalar_pdt[singular_indices[i-1]] - singular_scalar_pdt[singular_indices[i]]) <= std::numeric_limits<double>::epsilon()){
                        _singular_values[len] = _singular_values[len] + crit_mul;
                    }else{
                        _singular_T.push_back(singular_scalar_pdt[singular_indices[i]]);
                        len++;
                        _singular_values.push_back(_Values[find_euler_car_index(_T,_singular_T[len],0,_T.size())] - cla_crit_val[index][singular_indices[len]]);
                    }
                }
            }

            Radon_transform radon_transform(_T, _Values, singular_scalar_pdt, _singular_values);
            return radon_transform;
        }

        std::size_t find_euler_car_index(std::vector<double> &table, double t, std::size_t begin, std::size_t end){
            if(end - begin <= 1){
                return begin;
            }else{
                std::size_t index = (begin + end) / 2;
                if(table[index] > t){
                    return find_euler_car_index(table, t, begin, index);
                }else{
                    return find_euler_car_index(table, t, index, end);
                }
            }
            
        }

        std::vector<std::vector<double>> compute_radon_transform_python(std::vector<double> &direction){
            Radon_transform radon = compute_radon_transform(direction);
            std::vector<std::vector<double>> result;

            result.push_back(radon.T);
            result.push_back(radon.Values);
            result.push_back(radon.singular_T);
            result.push_back(radon.singular_values);

            return result;
        }

        //*********************************************//
        //Functions to compute euler caracteristic transform
        //*********************************************//
        std::vector<int> principal_direction(std::vector<double> &v){
            std::vector<int> pdv;
            for(std::size_t i=0; i<v.size(); i++){
                if(v[i] > 0){
                    pdv.push_back(1);
                }else{
                    pdv.push_back(-1);
                }
            }
            return pdv;
        }

        std::vector<std::vector<double>> compute_ect_python(std::vector<double> &direction){
            Euler_caracteristic_transform ect = compute_ect(direction);
            std::vector<std::vector<double>> res;

            res.push_back(ect.T);
            res.push_back(ect.transform_values);

            return res;
        }

        Euler_caracteristic_transform compute_ect(std::vector<double> &direction){
            // To adapt
            if(are_cla_crit_pts_computed == 0){
                compute_classical_critical_vertices();
            }

            int index = get_vector_index(direction);
            if(cla_crit_pts[index].size() == 0){
                std::vector<double> tmp;
                tmp.push_back(0);
                Euler_caracteristic_transform ect(tmp, tmp);
                return ect;
            }
            
            std::vector<double> scalar_pdt;
            std::vector<int> indices;

            for(std::size_t i=0; i<cla_crit_pts[index].size(); i++){
                scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[cla_crit_pts[index][i]]].begin(), 0.0));
                indices.push_back(i);
            }
            std::sort(indices.begin(), indices.end(), [&scalar_pdt](int i, int j) {return scalar_pdt[i] < scalar_pdt[j];});

            std::vector<double> sorted_scalar_products;
            std::vector<double> euler_car_accumulator;
            sorted_scalar_products.push_back(scalar_pdt[indices[0]]);
            euler_car_accumulator.push_back(cla_crit_val[index][indices[0]]);       
            for(int i=1; i<(int)indices.size(); i++){
                if(std::abs(scalar_pdt[indices[i-1]] - scalar_pdt[indices[i]]) <= std::numeric_limits<double>::epsilon()){
                    euler_car_accumulator[euler_car_accumulator.size()-1] += cla_crit_val[index][indices[i]];
                }else{
                    sorted_scalar_products.push_back(scalar_pdt[indices[i]]);
                    euler_car_accumulator.push_back(cla_crit_val[index][indices[i]] + euler_car_accumulator[euler_car_accumulator.size()-1]);
                }
            }
            // Passer par ref dans la classe Euler_car..
            Euler_caracteristic_transform ect(sorted_scalar_products, euler_car_accumulator);
            return ect;
        }

        int compute_sum_dimcell(){
            int s = 0;
            for(Simplex_key i=0; i<this->num_simplices(); i++){
                s += this->dimension(i)+1;
            }
            return s;
        }

        // Functions to get attributes
        std::vector<int> get_ordinary_critical_vertices(int index){
            return ord_crit_pts[index];
        }

        std::vector<int> get_ordinary_critical_values(int index){
            return ord_crit_val[index];
        }

        std::vector<int> get_classical_critical_vertices(int index){
            return cla_crit_pts[index];
        }

        std::vector<int> get_classical_critical_values(int index){
            return cla_crit_val[index];
        }

        std::vector<double> get_vertex_coordinates(int index){
            return embedding[embedding_index[index]];
        }

        // TODO DELETE : Warning - IN WRAPPER
        double compute_euler_caracteristic_of_complex(){
            double s = 0;
            for(Simplex_key i=0; i<this->num_simplices(); i++){
                if(this->dimension(i) % 2 == 0){
                    s += this->filtration(i);
                }else{
                    s -= this->filtration(i);
                }
                
            }
            return s;
        }

};

#endif // EMBEDDED_CUBICAL_COMPLEX_H_