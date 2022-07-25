### **BENCHMARKING**

The following program is able to test the computation time of hybrid transforms on your computer in two ways :
- Compare the different methods using a single thread
- Compare the computation time as a function of the number of threads used.

## **Compilation**
To compile the benchmarking program just use :
```
g++ -pthread benchmarking.cpp -o program_name
```

## **Usage**
Once you've compiled the program you can use it with the following arguments :

./program_name b n s<sub>0</sub> s<sub>1</sub> ... s<sub>n-1</sub> r p k m 

**Arguments :**
- b : if zero, runs all the methods on one core, if 1 runs the critical points method with 1,2,3... cores up to the number of cores that you have.
- n : the dimension of the complexes.
- s<sub>0</sub>,...,s<sub>n-1</sub> : the sizes of your complex.
- r : the complexes will have values in [|0,r_1|].
- p : the number of complexes to create. 
- k : the vectors will have their coefficients in [-k,k].
- m : the number of vectors to create for each complex.


The time displayed by the program is the real elapsed time.
