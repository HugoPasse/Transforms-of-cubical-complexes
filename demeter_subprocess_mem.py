import numpy as np
from test_shapes import regular_points

import demeter.euler as euler
import demeter.misc as misc
import demeter.directions as dirs

import sys
import time
import os, psutil
process = psutil.Process(os.getpid())

if len(sys.argv) != 7:
	print('Incorrect arguments')
	exit()

n = sys.argv[1] 
spacing = sys.argv[2]
n_dir = sys.argv[3]
dim = sys.argv[4]
path_to_savings = sys.argv[5]
T = sys.argv[6]


if not (n.isdigit() and spacing.isdigit() and dual.isdigit() and T.isdigit()):
	print("Arguments must be integers")
	exit()

if not path_to_savings.isprintable():
	print("path_to_savings must be printable")
	exit()

n = int(n)
spacing = int(spacing)
n_dir = int(n_dir)
dual = int(dual)
dim = int(dim)
T = int(T)
if not dual == 0 : dual = 1

img = np.array(regular_points((n,n),np.array([spacing,spacing]),np.array([spacing,spacing])))

# Initialize complex
mtic = process.memory_full_info().uss
tic = time.perf_counter()
cplx = euler.CubicalComplex(img)
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_cplx = toc-tic
print(f'{time_cplx} {(mtoc-mtic)/1000}', end=' ')

# Complexifying
mtic = process.memory_full_info().uss
tic = time.perf_counter()
seed = cplx.complexify()
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_cplxying = toc-tic
print(f'{time_cplxying} {(mtoc-mtic)/1000}', end=' ')

# Computing transform
directions = np.random.rand(n_dir, dim)
mtic = process.memory_full_info().uss
tic = time.perf_counter()
value = seed.ECT(directions, T)
toc = time.perf_counter()
mtoc = process.memory_full_info().uss
time_transform = toc-tic
print(f'{time_transform} {(mtoc-mtic)/1000}')