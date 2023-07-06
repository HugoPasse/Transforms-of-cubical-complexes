import numpy as np
import embedded_cubical_complex as ecc
from test_shapes import regular_points
import sys
import time
import os, psutil
process = psutil.Process()

if len(sys.argv) != 5:
	print('Incorrect arguments')
	exit()

n = sys.argv[1]
spacing = sys.argv[2]
dual = sys.argv[3]
path_to_savings = sys.argv[4]

if not (n.isdigit() and spacing.isdigit() and dual.isdigit()):
	print("Arguments must be integers")
	exit()

if not path_to_savings.isprintable():
	print("path_to_savings must be printable")
	exit()

n = int(n)
spacing = int(spacing)
dual = int(dual)
if not dual == 0 : dual = 1

img = np.array(regular_points((n,n),np.array([spacing,spacing]),np.array([spacing,spacing])))

# Construct complex
tic = time.perf_counter()
mtic = process.memory_info().rss
cplx = ecc.EmbeddedComplex(img*(-1)**dual, 1, input_top_cells=(not dual))
mtoc = process.memory_info().rss
toc = time.perf_counter()
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
mtic = process.memory_info().rss
tic = time.perf_counter()
cplx.init_ect(1)
toc = time.perf_counter()
mtoc = process.memory_info().rss
time_preproc = toc-tic
print(f'{time_preproc} {(mtoc-mtic)/1000}', end=' ')

# Computing transform
transform = 0
mtic = process.memory_info().rss
tic = time.perf_counter()
for i in range(10):
	transform = cplx.compute_euler_caracteristic_transform(np.random.rand(2))
toc = time.perf_counter()
mtoc = process.memory_info().rss
time_transform = toc-tic
print(f'{time_transform} {(mtoc-mtic)/1000}', end=' ')


cplx.init_hybrid_transform(1)
print(f'{len(cplx.get_critical_vertices(0))}')