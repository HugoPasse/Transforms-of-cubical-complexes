import numpy as np
import embedded_cubical_complex as ecc
from test_shapes import regular_points
import sys
import time
import os, psutil
process = psutil.Process(os.getpid())

if len(sys.argv) != 8:
	print('Incorrect arguments')
	exit()

n = sys.argv[1] 
spacing = sys.argv[2]
n_dir = sys.argv[3]
dual = sys.argv[4]
dim = sys.argv[5]
path_to_savings = sys.argv[6]
transform = sys.argv[7]


if not (n.isdigit() and spacing.isdigit() and dual.isdigit()):
	print("Arguments must be integers")
	exit()

if not path_to_savings.isprintable():
	print("path_to_savings must be printable")
	exit()

if not transform in ['ECT', 'Radon', 'HT']:
	print('transform must be either ECT, Radon or HT')
	exit()

n = int(n)
spacing = int(spacing)
n_dir = int(n_dir)
dual = int(dual)
dim = int(dim)
if not dual == 0 : dual = 1

img = np.array(regular_points((n,n),np.array([spacing,spacing]),np.array([spacing,spacing])))

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
directions = np.random.rand(n_dir, dim)
mtic = process.memory_full_info().uss
tic = time.perf_counter()
if transform == 'ECT':
	for i in range(n_dir):
		value = cplx.compute_euler_caracteristic_transform(np.random.rand(dim))
elif transform == 'Radon':
	for i in range(n_dir):
		value = cplx.compute_radon_transform(np.random.rand(dim))
elif transform == 'HT':
	value = cplx.compute_hybrid_transform('cos', directions, num_jobs=1)
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_transform = toc-tic
print(f'{time_transform} {(mtoc-mtic)/1000}', end=' ')

cplx.preproc_radon_transform()
print(f'{len(cplx.get_classical_critical_vertices(3))}', end=' ')
print(f'{len(cplx.get_ordinary_critical_vertices(3))}', end=' ')
