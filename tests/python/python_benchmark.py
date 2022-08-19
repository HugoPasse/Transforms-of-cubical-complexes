import embedded_cubical_complex as embedded_cc

import numpy as np
import time

def create_random_complex(sizes, val_range):
	return np.random.randint(val_range, size=sizes)

def random_vectors(number, dimension, val_range):
	return (np.multiply(np.random.rand(number,dimension), val_range))

def benchmark(iterations, num_cplx, cplx_sizes, vect_range, cplx_range, num_cores=0):
	dimension = len(cplx_sizes)
	total_time = 0
	vectors = random_vectors(iterations, dimension, vect_range)
	tot_pc = 0
	tot_c = 0
	for i in range(num_cplx):
		rdcplx = create_random_complex(cplx_sizes,cplx_range)
		prec = time.time()
		cplx = embedded_cc.EmbeddedComplex(rdcplx, num_cores)
		cplx.init_hybrid_transform(num_cores)
		start = time.time()
		T = cplx.compute_hybrid_transform("exp", vectors, num_cores	)
		end = time.time()
		tot_c += end - start
		tot_pc += start - prec
		print(T)
	print("Precalc time : " + str(tot_pc) + ", computation time : " + str(tot_c) + ", num_cores : " + str(num_cores))

for i in range(1,9):
	for j in range(10):
		benchmark(10000,1,np.array([100,100]),10,255,i)