#include <gudhi/Bitmap_cubical_complex.h>

#include <future>
#include <thread>
#include <chrono>

#ifndef EMBEDDED_CPLX
#define EMBEDDED_CPLX

template <typename T>
void print_vector(std::vector<T> vect){
    if(vect.size() == 0){
        std::cout << "[]\n";
    }else{
        std::cout << "[" << vect[0];
        for(long unsigned int i=1; i<vect.size(); i++){
            std::cout << ", " << vect[i];
        }
        std::cout << "]\n";
    }
}

template <typename T>
class Embedded_cubical_complex : public Gudhi::cubical_complex::Bitmap_cubical_complex<T>
{

    private:

        typedef std::size_t Simplex_key;
        typedef typename T::filtration_type Filtration_value;
        typedef Simplex_key Simplex_handle;

    public:
        std::vector<std::vector<double>> embedding;         //Array of the cubical complexe's points' coordinates
        std::vector<int> embedding_index;                   //Array to link vertices index in the cubical complex to their index in the embedding

        std::vector<Simplex_handle> relevant_index;         //An optimization for sparse images, we only compute the sum on non-zero indexes

        std::vector<std::vector<int>> min_arranging;        //Arranging vectors, given a vector e = (e_0, ... , e_d-1) you can know for each simplex given by its key k the index of its vertices that verify the min / max of the scalar product with e
        std::vector<std::vector<int>> max_arranging;        //min_index = min_arranging[index[key]] where index is the sum of 2^i with i such that e_i is positive

        std::vector<int> sizes_pdt;                         //Products of the sizes from s_0 up to s_0*s_1*...*s_(n-2)

        //********* Attributes for critical points optim *********//
        std::vector<std::vector<int>> critical_vertices_old;    //critical_vertices[index] gives you the vector of the critical points' indexes that have e for principal direction, where e is encoded by i with the following encoding :
                                                                //e = (e_1, ..., e_n-1), index = sum of 2^i with i such that e_i is positive.
                                                                //The vector of the critical points you have to substract to the sum is critical_vertices[2^dimension - 1 - index]
        //*********************************************//
        std::vector<std::vector<int>> critical_vertices;
        std::vector<std::vector<int>> critical_multiplicity;

        //*********************************************//
        //Constructors
        //*********************************************//
        Embedded_cubical_complex(const std::vector<unsigned>& dimensions,
            const std::vector<Filtration_value>& top_dimensional_cells):Gudhi::cubical_complex::Bitmap_cubical_complex<T>(dimensions, top_dimensional_cells)
            {
                sizes_pdt.push_back(2*this->sizes[0]+1);
                for(Simplex_handle i=1; i<this->dimension()-1; i++){           //In this loop we compute the product of the number of cubes in each direction, it optimizes the complexity of the get_coordinates_in_complex function
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
            for(Simplex_handle i = 0; i<this->num_simplices();i++){
                if(this->dimension(i) == 0){
                    std::vector<int> coords = get_coordinates_in_complex(i);
                    std::vector<double> embedded_coords(coords.size());
                    for(Simplex_handle j=0; j<coords.size(); j++){
                        embedded_coords[j] = (double)(coords[j] / 2) / this->sizes[j] - 0.5;
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

        void initalize_relevant_indexes(){
            for(Simplex_handle handle = 0; handle < this->num_simplices(); handle++){
                if(this->filtration(handle) != 0 && this->dimension(handle) != 0){
                    relevant_index.push_back(handle);
                }
            }
        }

        void compute_arrangement(){
            int dim = this->dimension();
            std::vector<int> e(dim,-1);

            for(int i=0; i<std::pow(2,dim); i++){       //Loop on the 2^d possible directions for e

                std::vector<int> min_vertices;          //Vectors to store, given e, the indexes of min and max vertex for each simplex
                std::vector<int> max_vertices;          
                for(Simplex_handle key = 0; key<this->num_simplices(); key++){     //Loop on the simplexes
                    std::vector<int> vertices = get_cell_vertices(key);

                    double min = __DBL_MAX__;
                    double max = -__DBL_MAX__;

                    int index_min = -1;
                    int index_max = -1;

                    for(int v=0; v<(int)vertices.size(); v++){           //Loop on the vertices of the simplex
                        double pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[vertices[v]]].begin(),0.0);
                        if(pdt < min){
                            min = pdt;
                            index_min = vertices[v];
                        }
                        if(pdt > max){
                            max = pdt;
                            index_max = vertices[v];
                        }
                    }

                    min_vertices.push_back(index_min);
                    max_vertices.push_back(index_max);
                }

                for(int j=0; j<dim; j++){       //Changing direction vector
                    if(e[j] == 1){
                        e[j] = -1;
                    }else{
                        e[j] = 1;
                        break;
                    }
                }
                min_arranging.push_back(min_vertices);
                max_arranging.push_back(max_vertices);
            }
        }

        //*********************************************//
        //Functions for critical points
        //*********************************************//

        void compute_non_singular_critical_vertices(int num_jobs = 0){
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
                
                if(multiplicity != 0){                                  //If relevant we add it to the critical points
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

        int compute_euler_car_in_direction(Simplex_handle vertex, std::vector<int> direction, int reverse_vector){ //Return euler car with multiplicity of the intersection between the cells in the neigbourhood of the vertex and the hyperplane orthogonal to direction in the neighbourhood of the vertex
            int euler_car = 0;
            
            std::vector<int> coordinates = get_coordinates_in_complex(vertex);
            
            int dim = direction.size();
            std::vector<int> tmp(dim);
            int simplex_dim_sign = 1;
            
            for(int i=0; i<std::pow(2,dim)-1; i++){  //Looping on all of the adjacent cell in direction
                int k = 0;
                while(tmp[k] == 1){
                    tmp[k] = 0;
                    coordinates[k] -= reverse_vector * direction[k];
                    simplex_dim_sign *= -1;
                    k++;
                }
                coordinates[k] += reverse_vector * direction[k];

                if(k < dim){
                    tmp[k] = 1;
                    simplex_dim_sign *= -1;
                }
                if(are_coordinates_in_complex(coordinates) == 1){
                    Simplex_key key = get_key_from_coordinates(coordinates);
                    euler_car += this->filtration(key) * simplex_dim_sign;
                }
            }

            return euler_car;
        }


        //*********************************************//
        //Printing functions
        //*********************************************//

        void print_embedding(){
            std::cout << "[";
            for(Simplex_handle i=0; i<embedding_index.size(); i++){
                if(embedding_index[i] != -1){
                    print_vector(embedding[embedding_index[i]]);
                }
            }
            std::cout << "]\n";
        }

        void print_arrangement(){
            std::cout << "MIN : [";
            for(Simplex_handle i=0; i<min_arranging.size(); i++){
                print_vector(min_arranging[i]);
            }
            std::cout << "]\n";
            std::cout << "MAX : [";
            for(Simplex_handle i=0; i<max_arranging.size(); i++){
                print_vector(max_arranging[i]);
            }
            std::cout << "]\n";
        }

        void print_critical_vertices_old(){
            std::cout << "Critical vertices : \n";
            for(std::size_t i=0; i<critical_vertices_old.size(); i++){
                print_vector(critical_vertices_old[(int)i]);
            }
            std::cout << "\n";
        }

        void print_critical_vertices(){
            std::cout << "Critical vertices test : \n";
            for(std::size_t i=0; i<critical_vertices.size(); i++){
                print_vector(critical_vertices[(int)i]);
            }
            std::cout << "\n";
        }

        void print_critical_multiplicity(){
            std::cout << "Critical multiplicity : \n";
            for(std::size_t i=0; i<critical_multiplicity.size(); i++){
                print_vector(critical_multiplicity[(int)i]);
            }
            std::cout << "\n";
        }

        void print_filtration(){
            std::cout << "Filtration : \n[";
            int n = this->num_simplices()-1;
            for(int i=2*this->sizes[1]; i>=0; i--){
                for(int j=0; j<2*this->sizes[0]+1; j++){
                    std::cout << this->filtration(j+i*(2*this->sizes[0]+1)) << ", ";
                }
                std::cout << "\n";
            }
            std::cout << "]\n";
        }

        //*********************************************//
        //Those are arithmetic functions to find the index of the vertex within the vertices
        //*********************************************//

        //This function gives the coordinates of the cell given by key
        std::vector<int> get_coordinates_in_complex(Simplex_key key){

            int n = (int)this->dimension();

            std::vector<int> coordinates;      //Coordinates of the face indexed by key to return
            
            for(int i=n-1; i>0; i--){
                //std::cout << "i = " << i << " key = " << key << " sizes_pdt[i] = " << sizes_pdt[i-1] << "\n";
                coordinates.insert(coordinates.begin(),key / sizes_pdt[i-1]);
                key = key - coordinates[0]*sizes_pdt[i-1];
            }

            coordinates.insert(coordinates.begin(),key);

            return coordinates;
        }

        //This functions gives you the key of the vertex given the coordinates of it
        Simplex_key get_key_from_coordinates(std::vector<int> coordinates){
            Simplex_key key = 0;
            for(int i=this->dimension()-1; i>=0; i--){
                key = key*(2*this->sizes[i] + 1) + coordinates[i];
            }
            return key;
        }

        int get_simplex_embedding_index(Simplex_key key){                       //THIS FUNCTION MIGHT BE USLESS AS IT IS SLOW COMPARED TO ACCESSING THE VECTOR "embedding_index"
            std::vector<int> coordinates = get_coordinates_in_complex(key);     //Get the coordinates in the complex

            Simplex_key embedding_key = 0;

            for(int i=this->dimension()-1; i>=0; i--){
                embedding_key = (this->sizes[i]+1)*embedding_key + (coordinates[i] / 2);
            }

            return embedding_key;
        }

        //This function returns a vector with the keys of the vertices of the cell given by key
        std::vector<int> get_cell_vertices(Simplex_key key){
            std::vector<int> cell_coordinates = get_coordinates_in_complex(key);    //Get the coordinates of cel indexed by key
            int n = cell_coordinates.size();

            std::vector<int> odd_coordinates;
            std::vector<int> cell_vertex(n);
            int n_odd_coords = 0;
            
            for(int i=0; i<n; i++){                         //computring the number and indexes of odd coordinates in vector cell_coordinates
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

            for(int i=0; i<std::pow(2,n_odd_coords)-1;i++){
                for(int j=0; j<n_odd_coords;j++){                   //We use a binary counter on the odd cordinates of the simplex to compute all of the vertices coordinates
                    if(tmp_vect[j] == 1){
                        tmp_vect[j] = 0;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] - 2;
                    }else{
                        tmp_vect[j] = 1;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] + 2;
                        break;
                    }
                }
                
                cell_vertices.push_back(get_key_from_coordinates(cell_vertex));     //Adding the key of the vertice we found in cell_vertices
            }

            return cell_vertices;
        }

        //Given a direction vector e, return the index in the arrangement corresponding to the vector of vertices that verifies the min and the max of the scalar product within each cell of the complex
        int get_vector_index(std::vector<double> e){
            int index = 0;
            int muliplier = 1;

            for(int i=0; i<(int)e.size(); i++){
                if(e[i] >= 0){
                    index += muliplier;
                }
                muliplier = muliplier*2;
            }

            return index;
        }

        int are_coordinates_in_complex(std::vector<int> coordinates){
            std::size_t coord_dim = coordinates.size();
            if(coord_dim != this->sizes.size()){
                return 0;
            }

            for(std::size_t i=0; i<coord_dim; i++){
                if(coordinates[i] < 0 || coordinates[i] > 2*((int)this->sizes[i])){
                    return 0;
                }
            }
            return 1;
        }

        //*********************************************//
        //Functions to compute hybrid transform, needs a functions to compute kernel's antiderivative (abusevely called kernel) and a direction vector e
        //*********************************************//

        double compute_hybrid_transform_simple(double (*kernel)(double), std::vector<double> e){

            std::vector<double> scalar_products(this->embedding.size(),0);       //Used to store scalar products
            std::vector<int> is_pdt_computed(this->embedding.size(),0);       //Used to know if scalar product has already been computed

            double sum = 0.0;

            for(Simplex_handle key=0; key<this->num_simplices(); key++){               //Looping on all of the cells
                if(this->dimension(key) > 0 && this->filtration(key) != 0){
                    std::vector<int> vertices = get_cell_vertices(key);
                    double min_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[vertices[0]]].begin(),0.0);
                    double max_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[vertices[0]]].begin(),0.0);

                    scalar_products[embedding_index[vertices[0]]] = min_pdt;
                    is_pdt_computed[embedding_index[vertices[0]]] = 1;

                    for(Simplex_handle v=1; v<vertices.size(); v++){   //Finding the min and max scalar products 
                        double pdt = 0;
                        if(is_pdt_computed[embedding_index[v]] == 1){
                            pdt = scalar_products[embedding_index[vertices[v]]];
                        }else{
                            pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[vertices[v]]].begin(),0.0);
                            
                            scalar_products[embedding_index[vertices[v]]] = pdt;
                            is_pdt_computed[embedding_index[vertices[v]]] = 1;
                        }
                        min_pdt = std::min(pdt,min_pdt);
                        max_pdt = std::max(pdt,max_pdt);
                    }

                    double sgn = 1.0;     
                    if(this->get_dimension_of_a_cell(key)%2 == 0){sgn = -1.0;}

                    sum += sgn * (double)this->filtration(key) * (kernel(max_pdt) - kernel(min_pdt));    //Adding a new term to the sum                
                }
            }
            return sum;
        }

        //*********************************************//
        // First optimized version with arrangements
        //*********************************************//
        double compute_hybrid_transform_arrangement(double (*kernel)(double), std::vector<double> e){
            std::vector<double> scalar_products(this->embedding.size(),0);       //Used to store scalar products
            std::vector<int> is_pdt_computed(this->embedding.size(),0);          //Used to know if scalar product has already been computed
            
            int index = get_vector_index(e);                                

            double sum = 0.0;

            double min_pdt = 0;
            double max_pdt = 0;

            for(Simplex_handle key=0; key<this->num_simplices(); key++){               //Looping on all of the cells
                if(this->dimension(key) > 0 && this->filtration(key) != 0){            
                    
                    if(is_pdt_computed[embedding_index[min_arranging[index][key]]] == 1){
                        min_pdt = scalar_products[embedding_index[min_arranging[index][key]]];
                    }else{
                        min_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[min_arranging[index][key]]].begin(),0.0);
                        scalar_products[embedding_index[min_arranging[index][key]]] = min_pdt;
                        is_pdt_computed[embedding_index[min_arranging[index][key]]] = 1;
                    }

                    if(is_pdt_computed[embedding_index[max_arranging[index][key]]] == 1){
                        max_pdt = scalar_products[embedding_index[max_arranging[index][key]]];
                    }else{
                        max_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[max_arranging[index][key]]].begin(),0.0);
                        scalar_products[embedding_index[max_arranging[index][key]]] = max_pdt;
                        is_pdt_computed[embedding_index[max_arranging[index][key]]] = 1;
                    }
                    
                    
                    double sgn = 1.0;     
                    if(this->get_dimension_of_a_cell(key)%2 == 0){sgn = -1.0;}

                    

                    sum += sgn * (double)this->filtration(key) * (kernel(max_pdt) - kernel(min_pdt));    //Adding a new term to the sum                
                }
            }
            return sum;
        }

        //*********************************************//
        // Optimized version arrangments only using relevant indexes
        //*********************************************//
        double compute_hybrid_transform_sparse(double (*kernel)(double), std::vector<double> e){
            std::vector<double> scalar_products(this->embedding.size(),0);       //Used to store scalar products
            std::vector<int> is_pdt_computed(this->embedding.size(),0);          //Used to know if scalar product has already been computed
            
            int index = get_vector_index(e);                                

            double sum = 0.0;

            for(Simplex_handle key=0; key<relevant_index.size(); key++){               //Looping only on relevant cells

                double min_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[min_arranging[index][relevant_index[key]]]].begin(),0.0); 
                double max_pdt = std::inner_product(e.begin(),e.end(),embedding[embedding_index[max_arranging[index][relevant_index[key]]]].begin(),0.0);
                
                double sgn = 1.0;     
                if(this->get_dimension_of_a_cell(relevant_index[key])%2 == 0){sgn = -1.0;}

                sum += sgn * (double)this->filtration(relevant_index[key]) * (kernel(max_pdt) - kernel(min_pdt));    //Adding a new term to the sum                
            
            }
            return sum;
        }

        //*********************************************//
        //  Using critical points to compute the transform  
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
};

#endif

/* CE CODE FONCTIONNE POUR LE MULTITHREADING

std::promise<double> promiseObj;
std::future<double> futureObj = promiseObj.get_future();

std::thread th(&Embedded_cubical_complex::launch_async_transform, this, &promiseObj, kernel, vect_list[0]);

while(futureObj.wait_for(wait_time) != std::future_status::ready){
    std::cout << "Waiting\n";
}

std::cout << "Result : " << futureObj.get() << "\n";
th.join();

*/