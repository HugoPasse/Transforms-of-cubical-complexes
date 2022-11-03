#include <gudhi/Bitmap_cubical_complex.h>

#include <future>
#include <thread>
#include <chrono>

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

        //Critical points for hybrid transforms (also used for radon transform)
        int are_non_singular_vertices_computed = 0;
        std::vector<std::vector<int>> critical_vertices;
        std::vector<std::vector<int>> critical_multiplicity;

        //Critical points for Radon transform
        int are_singular_vertices_computed = 0;
        std::vector<std::vector<int>> zero_measure_critical_vertices;
        std::vector<std::vector<int>> zero_measure_critical_multiplicity;

        //Points for Euler curve transform
        int are_ect_points_computed = 0;
        std::vector<std::vector<int>> ect_points;
        std::vector<std::vector<int>> ect_variations;

        //If we have a vector e = (e_1,...,e_n), its index is the sum of the 2^i were i are the indexes such that e_i >= 0
        //If index_v = index_w they have the same critical points. The critical points are stored in critical_vertices[index]

        //*********************************************//
        //Constructors
        //*********************************************//
        Embedded_cubical_complex(const std::vector<unsigned>& dimensions,
            const std::vector<Filtration_value>& top_dimensional_cells):Gudhi::cubical_complex::Bitmap_cubical_complex<T>(dimensions, top_dimensional_cells)
            {
                sizes_pdt.push_back(2*this->sizes[0]+1);
                //In this loop we compute the product of the number of cubes in each direction, it optimizes the complexity of the get_coordinates_in_complex function
                for(Simplex_handle i = 1; i < this->dimension()-1; i++){
                    sizes_pdt.push_back(sizes_pdt[i-1]*(2*this->sizes[i]+1));
                }
                initalize_embedding();
                initalize_embedding_index();
                impose_upper_star_filtration();
            }   

        void impose_upper_star_filtration(){
            for(Gudhi::cubical_complex::Bitmap_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>>::Top_dimensional_cells_iterator it = this->top_dimensional_cells_iterator_begin(); it != this->top_dimensional_cells_iterator_end(); ++it){
                std::vector<std::size_t> boundary = this->get_boundary_of_a_cell(*it);
                for(std::size_t i=0; i<boundary.size(); i++){
                    this->data[boundary[i]] = std::max(this->filtration(boundary[i]),this->filtration(*it));
                    impose_upper_star_filtration_from_simplex(boundary[i]);
                }
            }
        }

        void impose_upper_star_filtration_from_simplex(Simplex_handle sh){
            if(this->dimension(sh) > 0){
                std::vector<std::size_t> boundary = this->get_boundary_of_a_cell(sh);
                for(std::size_t i=0; i<boundary.size(); i++){
                    this->data[boundary[i]] = std::max(this->filtration(boundary[i]),this->filtration(sh));
                }
            }
        }

        //*********************************************//
        //Functions for pretratement
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

        void init_hybrid_transform(int num_jobs=0){
            compute_non_singular_critical_vertices(num_jobs);
        }

        void init_radon_transform(int num_jobs=0){
            compute_non_singular_critical_vertices(num_jobs);
            compute_singular_critical_vertices(num_jobs);
        }

        //*********************************************//
        //Functions for critical points
        //*********************************************//
        void compute_non_singular_critical_vertices(int num_jobs = 0){
                if(are_non_singular_vertices_computed == 0){
                if(num_jobs > (int)std::thread::hardware_concurrency() || num_jobs <= 0){
                    num_jobs = std::thread::hardware_concurrency();
                }
                num_jobs = std::min(num_jobs, (int)this->sizes[0]+1);    //Because at line 135, indices must designate vertices on the first line
                
                int dim = this->dimension();
                Simplex_key n = this->num_simplices();
                std::chrono::seconds zero_sec{0};

                int num_vertex = 1;
                for(std::size_t  i = 0; i < this->sizes.size(); i++){
                    num_vertex *= this->sizes[i] + 1;
                }

                std::vector<int> direction(dim,-1);
                long unsigned int index = 0;

                for(int j = 0; j < (1 << dim); j++){    //Loop on every possible critical direction 
                    std::vector<int> tmp;
                    critical_vertices.push_back(tmp);
                    critical_multiplicity.push_back(tmp);

                    std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                    std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                    std::vector<std::thread> thread_vector;

                    for(int i = 0; i < num_jobs; i++){
                        std::promise<int> promiseObj;
                        promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                        future_vect.push_back(promise_vect[i].get_future());
                        //Objects to get return values from the thread
                        thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_non_singular_critical_vertices_subroutine, this, std::move(promise_vect[i]), direction, dim, n, i, num_jobs));
                        thread_vector[i].detach();  //Thread is now running concurently
                    }

                    int job = 0;
                    while(job < num_jobs){
                        //Getting thread response and appending it to the corresponding vectors
                        if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                            std::vector<std::vector<int>> thread_res = future_vect[job].get();

                            critical_vertices[index].insert(critical_vertices[index].end(), thread_res[0].begin(), thread_res[0].end());
                            critical_multiplicity[index].insert(critical_multiplicity[index].end(), thread_res[1].begin(), thread_res[1].end());
                            job++;
                        }
                    }
                    //Finding next direction vector
                    int k=0;
                    while(direction[k] == 1 && k < dim){
                        direction[k] = -1;
                        index -= (1u << k);
                        k++;
                    }
                    direction[k] = 1;
                    index += (1u << k);
                }
            }
            are_non_singular_vertices_computed = 1;
        }

        void compute_non_singular_critical_vertices_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> crit_subroutine;
            std::vector<int> mult_subroutine;
            
            std::vector<int> coords(dim,0);     //On each pass throught the loop, it will store the coordinates of the simplex
            coords[0] = 2*job_index;
            Simplex_key vertex = get_key_from_coordinates(coords);   //Will store the index of the vertex

            while(vertex < n_simplices){  //Loop on all vertices
                
                int euler_1 = compute_euler_car_in_direction(vertex,direction,1);
                int euler_2 = compute_euler_car_in_direction(vertex,direction,-1);
                int multiplicity = euler_1 - euler_2;    //Computing critical point multiplicity in direction
                
                if(multiplicity != 0){                                      //If relevant we add it to the critical points
                    crit_subroutine.push_back(vertex);     
                    mult_subroutine.push_back(multiplicity);
                }

                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if((unsigned)coords[j] > 2*this->sizes[j]+1){
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
            res.push_back(crit_subroutine);
            res.push_back(mult_subroutine);
            promiseObj.set_value(res);
        }

        void compute_singular_critical_vertices(int num_jobs=0){
            if(are_singular_vertices_computed == 0){
                if(num_jobs > (int)std::thread::hardware_concurrency() || num_jobs <= 0){
                    num_jobs = std::thread::hardware_concurrency();
                }
                num_jobs = std::min(num_jobs, (int)this->sizes[0]+1);
                
                int dim = this->dimension();
                Simplex_key n = this->num_simplices();
                std::chrono::seconds zero_sec{0};

                int num_vertex = 1;
                for(std::size_t  i = 0; i < this->sizes.size(); i++){
                    num_vertex *= this->sizes[i] + 1;
                }

                std::vector<int> direction(dim,-1);
                long unsigned int index = 0;

                for(int j = 0; j < (1 << dim); j++){    //Loop on every possible critical direction 
                    std::vector<int> tmp;
                    zero_measure_critical_vertices.push_back(tmp);
                    zero_measure_critical_multiplicity.push_back(tmp);

                    std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                    std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                    std::vector<std::thread> thread_vector;

                    for(int i = 0; i < num_jobs; i++){
                        std::promise<int> promiseObj;
                        promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                        future_vect.push_back(promise_vect[i].get_future());
                        //Objects to get return values from the thread
                        thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_singular_critical_vertices_subroutine, this, std::move(promise_vect[i]), direction, dim, n, i, num_jobs));
                        thread_vector[i].detach();  //Thread is now running concurently
                    }

                    int job = 0;
                    while(job < num_jobs){
                        //Getting thread response and appending it to the corresponding vectors
                        if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                            std::vector<std::vector<int>> thread_res = future_vect[job].get();

                            zero_measure_critical_vertices[index].insert(zero_measure_critical_vertices[index].end(), thread_res[0].begin(), thread_res[0].end());
                            zero_measure_critical_multiplicity[index].insert(zero_measure_critical_multiplicity[index].end(), thread_res[1].begin(), thread_res[1].end());
                            job++;
                        }
                    }
                    //Finding next direction vector
                    int k=0;
                    while(direction[k] == 1 && k < dim){
                        direction[k] = -1;
                        index -= (1u << k);
                        k++;
                    }
                    direction[k] = 1;
                    index += (1u << k);
                }
            }
            are_singular_vertices_computed = 1;
        }

        void compute_singular_critical_vertices_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> crit_subroutine;
            std::vector<int> mult_subroutine;

            std::vector<int> coords(dim,0);     //On each pass throught the loop, it will store the coordinates of the simplex
            coords[0] = 2*job_index;
            Simplex_key vertex = get_key_from_coordinates(coords);   //Will store the index of the vertex

            while(vertex < n_simplices){  //Loop on all vertices
                
                int euler_1 = compute_euler_car_in_direction(vertex,direction,1);
                int euler_2 = compute_euler_car_in_direction(vertex,direction,-1);
                
               if(euler_1 == euler_2){
                    if(this->filtration(vertex) != -euler_1){
                        crit_subroutine.push_back(vertex);
                        mult_subroutine.push_back(-this->filtration(vertex) - euler_1);
                    }
                }else if(this->filtration(vertex) != -euler_1 && this->filtration(vertex) != -euler_2){
                    crit_subroutine.push_back(vertex);
                    mult_subroutine.push_back(-this->filtration(vertex) - euler_2);
                }

                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if((unsigned)coords[j] > 2*this->sizes[j]+1){
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
            res.push_back(crit_subroutine);
            res.push_back(mult_subroutine);
            promiseObj.set_value(res);
        }
        
        //Return euler car with multiplicity of the intersection between the cells in the neigbourhood of the vertex and the hyperplane orthogonal to direction in the neighbourhood of the vertex
        int compute_euler_car_in_direction(Simplex_handle vertex, std::vector<int> direction, int reverse_vector){
            int euler_car = 0;
            int dim = direction.size();

            std::vector<int> coordinates = get_coordinates_in_complex(vertex); //This vector will successively take the coordinates of adjacent cells involved in the Euler's caracteristic calculation
            
            std::vector<int> tmp(dim);  //This vector will help us to find all the adjacent cells involved in calculations
            int simplex_dim_sign = 1;
            
            for(int i = 0; i < (1 << dim)-1; i++){  //Looping on all of the adjacent cell in direction
                int k = 0;
                while(tmp[k] == 1){             //Finding the next adjacent cell 
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

        void init_ect(int num_jobs=0){
            if(are_ect_points_computed == 0){
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
                    ect_points.push_back(tmp);
                    ect_variations.push_back(tmp);

                    std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                    std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                    std::vector<std::thread> thread_vector;

                    for(int j = 0; j < num_jobs; j++){
                        std::promise<int> promiseObj;
                        promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                        future_vect.push_back(promise_vect[j].get_future());
                        //Objects to get return values from the thread
                        thread_vector.push_back(std::thread(&Embedded_cubical_complex::init_ect_subroutine, this, std::move(promise_vect[j]), direction, dimension, n_simplices, j, num_jobs));
                        thread_vector[j].detach();  //Thread is now running concurently
                    }
                    int job = 0;
                    std::chrono::seconds zero_sec{0};
                    while(job < num_jobs){
                        //Getting thread response and appending it to the corresponding vectors
                        if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                            std::vector<std::vector<int>> thread_res = future_vect[job].get();
                            ect_points[index].insert(ect_points[index].end(), thread_res[0].begin(), thread_res[0].end());
                            ect_variations[index].insert(ect_variations[index].end(), thread_res[1].begin(), thread_res[1].end());
                            job++;
                        }
                    }

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
            are_ect_points_computed = 1;
        }

        void init_ect_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> coords(dim,0);
            Simplex_key vertex = 2*job_index;

            coords[0] = 2*job_index;

            std::vector<int> subroutine_points;
            std::vector<int> subroutine_variations;

            while(vertex < n_simplices){        //Loop on all vertices
                int euler_car = compute_euler_car_in_direction(vertex, direction, -1) + this->filtration(vertex);
                if(euler_car != 0){
                    subroutine_points.push_back(vertex);
                    subroutine_variations.push_back(euler_car);
                }
                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if((unsigned)coords[j] > 2*this->sizes[j]+1){
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

        void print_critical_vertices(){
            std::cout << "Critical vertices : \n";
            for(std::size_t i = 0; i < critical_vertices.size(); i++){
                print_vector(critical_vertices[(int)i]);
            }
            std::cout << "\n";
        }

        void print_critical_multiplicity(){
            std::cout << "Critical multiplicity : \n";
            for(std::size_t i = 0; i < critical_multiplicity.size(); i++){
                print_vector(critical_multiplicity[(int)i]);
            }
            std::cout << "\n";
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

        int num_vertices(){
            int num = 1;
            for(int i=0; i<this->sizes.size(); i++){
                num *= this->sizes[i]+1;
            }
            return num;
        }

        //*********************************************//
        //Functions to compute hybrid transform
        //*********************************************//
        
        //Compute hybrid transform of the complex in direction e, kernel is the antiderivative of the real kernel of the transform
        double compute_hybrid_transform(double (*kernel)(double), std::vector<double> e){
            int index = get_vector_index(e);
            int reverse_vector = 1;
            //As multiplicity values are stored assuming that the last coordinate of direction is > 0, potentialy change index and set reverse_vector to -1
            if(index >= (int)critical_vertices.size()){
                reverse_vector = -1;
                index = (1 << e.size()) - 1 - index;
            }
            
            double sum = 0.0;

            for(std::size_t i = 0; i < critical_vertices[index].size(); i++){   //Looping on critical vertices
                sum +=  critical_multiplicity[index][i] * kernel(std::inner_product(e.begin(),e.end(),embedding[embedding_index[critical_vertices[index][i]]].begin(),0.0));
            }
            return reverse_vector * sum;
        }

        //An overload of previous function to support multithreading
        std::vector<double> compute_hybrid_transform(double (*kernel)(double), std::vector<std::vector<double>> vect_list, unsigned int num_threads = -1){
            if(are_non_singular_vertices_computed == 0){
                compute_non_singular_critical_vertices();
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
                
                thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_hybrid_transform_subvector, this, std::move(promise_vect[i]), kernel, vect_list, begin, end));
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
        std::vector<double> compute_hybrid_transform(int kernel_number, std::vector<std::vector<double>> vect_list, unsigned int num_threads = -1){
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
        void compute_hybrid_transform_subvector(std::promise<std::vector<double>> promiseObj, double (*kernel)(double), std::vector<std::vector<double>> vect_list, std::size_t begin_index, std::size_t end_index){
            std::vector<double> results;
            for(std::size_t i = begin_index; i < end_index; i++){
                results.push_back(compute_hybrid_transform(kernel,vect_list[i]));
            }
            promiseObj.set_value(results);
        }

        //*********************************************//
        //Functions to compute radon transform
        //*********************************************//
        
        Radon_transform compute_radon_transform(std::vector<double> direction){
            if(are_singular_vertices_computed == 0){
                compute_singular_critical_vertices();
            }
            if(are_non_singular_vertices_computed == 0){
                compute_non_singular_critical_vertices();
            }

            std::vector<double> _T;
            std::vector<double> _Values;

            std::vector<double> _singular_T;
            std::vector<double> _singular_values;

            int index = get_vector_index(direction);
            int reverse_multiplicity = -1;

            if(index >= (int)critical_vertices.size()){
                reverse_multiplicity = -1;
                index = (1 << direction.size()) - 1 - index;
            }

            std::vector<double> scalar_pdt;
            std::vector<double> singular_scalar_pdt;

            std::vector<std::size_t> indices;
            std::vector<std::size_t> singular_indices;
            
            //Computing scalar products with non singular critical vertices
            for(std::size_t i = 0; i < critical_vertices[index].size(); i++){
                scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[critical_vertices[index][i]]].begin(), 0.0));
                indices.push_back(i);
            }
            //Computing scalar products with singular critical vertices
            for(std::size_t i = 0; i < zero_measure_critical_vertices[index].size(); i++){
                singular_scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[zero_measure_critical_vertices[index][i]]].begin(), 0.0));
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
                euler_car = reverse_multiplicity * critical_multiplicity[index][indices[0]];
                _T.push_back(scalar_pdt[indices[0]]);
                _Values.push_back(euler_car);
                for(std::size_t i = 1; i < indices.size(); i++){        
                    int crit_mul = reverse_multiplicity * critical_multiplicity[index][indices[i]];
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
                euler_car = _Values[find_euler_car_index(_T,_singular_T[0],0,_T.size())];
                _singular_values.push_back(euler_car - zero_measure_critical_multiplicity[index][singular_indices[0]]);
                for(std::size_t i = 1; i < singular_indices.size(); i++){
                    int crit_mul = zero_measure_critical_multiplicity[index][singular_indices[i]];
                    if(std::abs(singular_scalar_pdt[singular_indices[i-1]] - singular_scalar_pdt[singular_indices[i]]) <= std::numeric_limits<double>::epsilon()){
                        _singular_values[len] = _singular_values[len] + crit_mul;
                    }else{
                        _singular_T.push_back(singular_scalar_pdt[singular_indices[i]]);
                        len++;
                        _singular_values.push_back(_Values[find_euler_car_index(_T,_singular_T[len],0,_T.size())] - zero_measure_critical_multiplicity[index][singular_indices[len]]);
                    }
                }
            }
            
            Radon_transform radon_transform(_T, _Values, _singular_T, _singular_values);
            return radon_transform;
        }

        std::size_t find_euler_car_index(std::vector<double> table, double t, std::size_t begin, std::size_t end){
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

        std::vector<std::vector<double>> compute_radon_transform_python(std::vector<double> direction){
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
        std::vector<int> principal_direction(std::vector<double> v){
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

        /*
        Euler_caracteristic_transform compute_ect(std::vector<double> direction){
            std::vector<int> principal_direction_vector = principal_direction(direction);

            std::vector<double> scalar_pdt;
            std::vector<int> indices;

            int k = 0;
            for(std::size_t i=0; i<embedding_index.size(); i++){    //Computing scalar products
                if(embedding_index[i] >= 0){
                    scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[i]].begin(), 0.0));
                    indices.push_back(k);
                    k++;
                }
            }
            std::sort(indices.begin(), indices.end(), [&scalar_pdt](int i, int j) {return scalar_pdt[i] < scalar_pdt[j];});

            int vertex = 0;
            std::vector<int> coords(this->dimension());
            int n_simplices = this->num_simplices();

            std::vector<double> euler_car_of_vertices;
            while(vertex < n_simplices){        //Loop on all vertices
                euler_car_of_vertices.push_back(compute_euler_car_in_direction(vertex, principal_direction_vector, -1) + this->filtration(vertex));
             
                coords[0] = coords[0] + 2;
                for(std::size_t j = 0; j < this->dimension()-1; j++){
                    if((unsigned)coords[j] > 2*this->sizes[j]+1){
                        coords[j] = 0;
                        coords[j+1] = coords[j+1] + 2;
                    }else{
                        break;
                    }
                }
                vertex = get_key_from_coordinates(coords);
            }

            std::vector<double> sorted_scalar_products;
            std::vector<double> euler_car_accumulator;
            sorted_scalar_products.push_back(scalar_pdt[indices[0]]);
            euler_car_accumulator.push_back(euler_car_of_vertices[indices[0]]);       

            for(int i=1; i<(int)indices.size(); i++){
                if(std::abs(scalar_pdt[indices[i-1]] - scalar_pdt[indices[i]]) <= std::numeric_limits<double>::epsilon()){
                    euler_car_accumulator[euler_car_accumulator.size()-1] += euler_car_of_vertices[indices[i]];
                }else{
                    sorted_scalar_products.push_back(scalar_pdt[indices[i]]);
                    euler_car_accumulator.push_back(euler_car_of_vertices[indices[i]] + euler_car_accumulator[euler_car_accumulator.size()-1]);
                }
            }

            std::vector<double> reduced_sorted_scalar_products;
            std::vector<double> reduced_euler_car_accumulator;

            //Reducing the size of the vectors.
            if(sorted_scalar_products.size() > 0){
                reduced_sorted_scalar_products.push_back(sorted_scalar_products[0]);
                reduced_euler_car_accumulator.push_back(euler_car_accumulator[0]);

                for(int i=1; i<(int)euler_car_accumulator.size(); i++){
                    if(std::abs(euler_car_accumulator[i] - euler_car_accumulator[i-1]) > std::numeric_limits<double>::epsilon()){
                        reduced_euler_car_accumulator.push_back(euler_car_accumulator[i]);
                        reduced_sorted_scalar_products.push_back(sorted_scalar_products[i]);
                    }
                }
            }
            
            Euler_caracteristic_transform ect(reduced_sorted_scalar_products, reduced_euler_car_accumulator);
            return ect;
        }*/

        std::vector<std::vector<double>> compute_ect_python(std::vector<double> direction){
            Euler_caracteristic_transform ect = compute_ect(direction);
            std::vector<std::vector<double>> res;

            res.push_back(ect.T);
            res.push_back(ect.transform_values);

            return res;
        }

        Euler_caracteristic_transform compute_ect(std::vector<double> direction){
            if(are_ect_points_computed == 0){
                init_ect();
            }

            int index = get_vector_index(direction);

            std::vector<double> scalar_pdt;
            std::vector<int> indices;

            for(std::size_t i=0; i<ect_points[index].size(); i++){
                scalar_pdt.push_back(std::inner_product(direction.begin(), direction.end(), embedding[embedding_index[ect_points[index][i]]].begin(), 0.0));
                indices.push_back(i);
            }
            std::sort(indices.begin(), indices.end(), [&scalar_pdt](int i, int j) {return scalar_pdt[i] < scalar_pdt[j];});

            std::vector<double> sorted_scalar_products;
            std::vector<double> euler_car_accumulator;
            sorted_scalar_products.push_back(scalar_pdt[indices[0]]);
            euler_car_accumulator.push_back(ect_variations[index][indices[0]]);       

            for(int i=1; i<(int)indices.size(); i++){
                if(std::abs(scalar_pdt[indices[i-1]] - scalar_pdt[indices[i]]) <= std::numeric_limits<double>::epsilon()){
                    euler_car_accumulator[euler_car_accumulator.size()-1] += ect_variations[index][indices[i]];
                }else{
                    sorted_scalar_products.push_back(scalar_pdt[indices[i]]);
                    euler_car_accumulator.push_back(ect_variations[index][indices[i]] + euler_car_accumulator[euler_car_accumulator.size()-1]);
                }
            }
            Euler_caracteristic_transform ect(sorted_scalar_products, euler_car_accumulator);
            return ect;
        }

        //TMP function
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

        // Functions to get attributes
        std::vector<int> get_critical_vertices(int index){
            return critical_vertices[index];
        }

        std::vector<int> get_critical_multiplicity(int index){
            return critical_multiplicity[index];
        }

        std::vector<double> get_vertex_coordinates(int index){
            return embedding[embedding_index[index]];
        }
};

#endif // EMBEDDED_CUBICAL_COMPLEX_H_