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
T = sys.argv[5]
full = sys.argv[6]
index = sys.argv[7]

n = int(n)
n_dir = int(n_dir)
dual = int(dual)
full = int(full)
T = int(T)
index = int(index) + 1
if not dual == 0 : dual = 1
if not full == 0 : full = 1

data = pd.read_csv('fashion_mnist/fashion-mnist_train.csv', skiprows=lambda x:x not in [0,index])

time_cplx, time_preproc, time_transform, mem_cplx, mem_preproc, mem_transform = 0, 0, 0, 0, 0, 0

if full:
	img = (data.iloc[:,1:]).to_numpy().reshape((28,28))
	values = np.unique(img)
	for val_pix in values:
		img = np.zeros_like(X)
		img[X >  val_pix] = 1
		img[X <= val_pix] = 0
		# Construct complex
		mtic = process.memory_full_info().uss
		tic = time.perf_counter()
		cplx = euler.CubicalComplex(img)
		toc = time.perf_counter()
		mtoc = process.memory_full_info().uss
		time_cplx += toc-tic
		mem_cplx += mtoc-mtic

		# Preprocessing
		mtic = process.memory_full_info().uss
		tic = time.perf_counter()
		seed = cplx.complexify()
		toc = time.perf_counter()
		mtoc = process.memory_full_info().uss
		time_preproc += toc-tic
		mem_preproc += mtoc-mtic

		# Computing transform
		value = 0
		directions = np.random.rand(n_dir, 2)
		mtic = process.memory_full_info().uss
		tic = time.perf_counter()
		value = seed.ECT(directions, T)
		toc = time.perf_counter()
		mtoc = process.memory_full_info().uss
		time_transform += toc-tic
		mem_transform += mtoc-mtic
else:
	X = (data.iloc[:,1:]).to_numpy().reshape((28,28))
	values = np.unique(X)
	med = np.quantile(values, 0.5)
	img = np.zeros_like(X)
	img[X>med] = 1
	# Construct complex
	mtic = process.memory_full_info().uss
	tic = time.perf_counter()
	cplx = euler.CubicalComplex(img)
	toc = time.perf_counter()
	mtoc = process.memory_full_info().uss
	time_cplx += toc-tic
	mem_cplx += mtoc-mtic

	# Preprocessing
	mtic = process.memory_full_info().uss
	tic = time.perf_counter()
	seed = cplx.complexify()
	toc = time.perf_counter()
	mtoc = process.memory_full_info().uss
	time_preproc += toc-tic
	mem_preproc += mtoc-mtic

	# Computing transform
	value = 0
	directions = np.random.rand(n_dir, 2)
	mtic = process.memory_full_info().uss
	tic = time.perf_counter()
	value = seed.ECT(directions, T)
	toc = time.perf_counter()
	mtoc = process.memory_full_info().uss
	time_transform += toc-tic
	mem_transform += mtoc-mtic

print(f'{time_cplx} {mem_cplx/1000}', end=' ')
print(f'{time_preproc} {mem_preproc/1000}', end=' ')
print(f'{time_transform} {mem_transform/1000}')

