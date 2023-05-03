import embedded_cubical_complex as ecc
import test_shapes
import time

import matplotlib.pyplot as plt
import matplotlib.collections  as mc

import numpy as np

# Returns the parameters to give to test_shapes.regular_points for a certain 
# number of critical points given the size of the image and its dimension (not exact value)
def num_crit_to_spacing(size,dim,ncrit):
	p = (size/2)*((1/float(ncrit))**(1/float(dim)))
	return int(p)

# Varies the number of critical points in 2D of a n*n array
def timings_by_critical_points(n,pts,n_samples,n_dirs):
	T = np.zeros((len(pts),n_samples,3))
	N = np.zeros(len(pts))
	for i in range(len(pts)):
		p = pts[i]
		for j in range(n_samples):
			directions = np.random.rand(n_dirs,2)
			data = test_shapes.regular_points((n,n),np.array([p,p]),np.array([p,p]))
			# Construct complex
			tic = time.perf_counter()
			cplx = ecc.EmbeddedComplex(data,1)
			toc = time.perf_counter()
			T[i,j,0] = toc-tic
			# pre-processing
			tic = time.perf_counter()
			cplx.init_ect(1)
			toc = time.perf_counter()
			T[i,j,1] = toc-tic
			# Computing 1 transform
			euler_tr = 0
			tic = time.perf_counter()
			for d in directions:
				euler_tr = cplx.compute_euler_caracteristic_transform(d)
			toc = time.perf_counter()
			T[i,j,2] = (toc-tic)/float(n_dirs)
			N[i] = len(cplx.get_ect_points(0))
	return T,N

def plot_timings_by_critical_points(n,pts,n_samples,n_dirs):
	T,N = timings_by_critical_points(n,pts,n_samples,n_dirs)
	MT = np.mean(T,axis=1)
	fig, axs = plt.subplots(2, 2)
	axs[0, 0].plot(N,MT[:,0])
	axs[0, 0].set_title('Complex creation')
	axs[0, 1].plot(N,MT[:,1])
	axs[0, 1].set_title('Pre-processing')
	axs[1, 0].plot(N,MT[:,2])
	axs[1, 0].set_title('Transform computation')
	plt.show()


# Constant number of critical points in a variable size image
def timings_by_size(sizes,spacings,n_samples,n_dirs):
	T = np.zeros((len(sizes),n_samples,3))
	C = np.zeros(len(sizes))
	N = np.zeros(len(sizes))
	for i in range(len(sizes)):
		p = spacings[i]
		n = sizes[i]
		for j in range(n_samples):
			directions = np.random.rand(n_dirs,2)
			data = test_shapes.regular_points((n,n),np.array([p,p]),np.array([p,p]))
			# Construct complex
			tic = time.perf_counter()
			cplx = ecc.EmbeddedComplex(data,1)
			toc = time.perf_counter()
			T[i,j,0] = toc-tic
			# pre-processing
			tic = time.perf_counter()
			cplx.init_ect(1)
			toc = time.perf_counter()
			T[i,j,1] = toc-tic
			# Computing 1 transform
			euler_tr = 0
			tic = time.perf_counter()
			for d in directions:
				euler_tr = cplx.compute_euler_caracteristic_transform(d)
			toc = time.perf_counter()
			T[i,j,2] = (toc-tic)/float(n_dirs)
			C[i] = len(cplx.get_ect_points(0)) # To Change 
			N[i] = cplx.compute_sum_dimcell()
	return T,C,N
	
def plot_timings_by_size(sizes,spacings,n_samples,n_dirs):
	T,C,N = timings_by_size(sizes,spacings,n_samples,n_dirs)
	print(N)
	MT = np.mean(T,axis=1)
	fig, axs = plt.subplots(2, 2)
	axs[0, 0].plot(sizes,MT[:,0])
	axs[0, 0].set_title('Complex creation')
	axs[0, 1].plot(sizes,MT[:,1])
	axs[0, 1].set_title('Pre-processing')
	X,YCC,YPC = [],[],[]
	for i in range(len(T[:,:,1])):
		for j in range(len(T[:,:,1][i])):
			YCC.append(T[:,:,0][i][j])
			YPC.append(T[:,:,1][i][j])
			X.append(N[i])
	axs[1, 0].scatter(X,YCC)
	axs[1, 0].set_title('Complex creation time over complexity')
	axs[1, 1].scatter(X,YPC)
	axs[1, 1].set_title('Pre-processing time over complexity')
	plt.show()


def timings_by_dimension(dim,size,spacing,n_samples,n_dirs):
	T = np.zeros(dim,n_samples,3)
	for i in range(dim):
		for j in range(n_samples):
			directions = np.random.rand(n_dirs,2)
			data = test_shapes.regular_points((size,size),np.array([spacing,spacing]),np.array([spacing,spacing]))
			# Construct complex
			tic = time.perf_counter()
			cplx = ecc.EmbeddedComplex(data,1)
			toc = time.perf_counter()
			T[i,j,0] = toc-tic
			# pre-processing
			tic = time.perf_counter()
			cplx.init_ect(1)
			toc = time.perf_counter()
			T[i,j,1] = toc-tic
			# Computing 1 transform
			tic = time.perf_counter()
			for d in directions:
				euler_tr = cplx.compute_euler_caracteristic_transform(d)
			toc = time.perf_counter()
			T[i,j,2] = (toc-tic)/float(n_dirs)
	return T,C,N

# Main computations
# plot_timings_by_critical_points(10000,[num_crit_to_spacing(10000,2,1000000*i) for i in range(1,11)],10,1000)
# S = [100*i for i in range(1,5)]
# plot_timings_by_size(S,[2**i for i in range(1,5)],10,100)