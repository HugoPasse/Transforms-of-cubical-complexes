import numpy as np
import embedded_cubical_complex as ecc
import sys
import time
import os, psutil
import pandas as pd
process = psutil.Process(os.getpid())

n = sys.argv[1] 
n_dir = sys.argv[2]
dual = sys.argv[3]
path_to_savings = sys.argv[4]
transform = sys.argv[5]
range_val = sys.argv[6]
index = sys.argv[7]

n = int(n)
n_dir = int(n_dir)
dual = int(dual)
num_thresholds = int(range_val)
index = int(index) + 1
if not dual == 0 : dual = 1

data = pd.read_csv('fashion_mnist/fashion-mnist_train.csv', skiprows=lambda x:x not in [0,index])
img = (data.iloc[:,1:]).to_numpy().reshape((28,28))

# X = (data.iloc[:,1:]).to_numpy().reshape((28,28))
# quantiles = [k/(num_thresholds+1) for k in range(1,num_thresholds+1)]
# values = np.unique(X)[1:]
# val_quantiles = [np.quantile(values, q) for q in quantiles]

# img = np.zeros_like(X)
# for val_pix in val_quantiles:
# 	img[X>val_pix] = val_pix

# Construct complex
mtic = process.memory_full_info().uss
tic = time.perf_counter()
cplx = ecc.EmbeddedComplex(img*(-1)**dual, 1, input_top_cells=(not dual))
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_cplx = toc-tic
print(f'{time_cplx} {(mtoc-mtic)/1000}', end=' ')

# Impose upper star
time_ustar = 0
if not dual:
	tic = time.perf_counter()
	cplx.impose_upper_star_filtration()
	toc = time.perf_counter()
	time_ustar = toc-tic
print(f'{time_ustar} 0', end=' ')

# Preprocessing
mtic = process.memory_full_info().uss
tic = time.perf_counter()
if transform == 'ECT':
	cplx.preproc_ect(1)
elif transform == 'Radon':
	cplx.preproc_radon_transform(1)
elif transform == 'HT':
	cplx.preproc_hybrid_transform(1)
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_preproc = toc-tic
print(f'{time_preproc} {(mtoc-mtic)/1000}', end=' ')

# Computing transform
value = 0
directions = np.random.rand(n_dir, 2)
mtic = process.memory_full_info().uss
tic = time.perf_counter()
if transform == 'ECT':
	for i in range(n_dir):
		value = cplx.compute_euler_caracteristic_transform(np.random.rand(2))
elif transform == 'Radon':
	for i in range(n_dir):
		value = cplx.compute_radon_transform(np.random.rand(2))
elif transform == 'HT':
	value = cplx.compute_hybrid_transform('cos', directions, num_jobs=1)
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_transform = toc-tic
print(f'{time_transform} {(mtoc-mtic)/1000}', end=' ')

cplx.preproc_radon_transform()
print(f'{len(cplx.get_classical_critical_vertices(3))}', end=' ')
print(f'{len(cplx.get_ordinary_critical_vertices(3))}', end=' ')
