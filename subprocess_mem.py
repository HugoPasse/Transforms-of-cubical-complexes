import numpy as np
import embedded_cubical_complex as ecc
from test_shapes import regular_points
import sys
import os, psutil
process = psutil.Process()

if len(sys.argv) != 4:
	print('Incorrect arguments')
	exit()

n = sys.argv[1]
spacing = sys.argv[2]
dual = sys.argv[3]

if not (n.isdigit() and spacing.isdigit() and dual.isdigit()):
	print("Arguments must be integers")
	exit()

n = int(n)
spacing = int(spacing)
dual = int(dual)
if not dual == 0 : dual = 1

img = np.array(regular_points((n,n),np.array([spacing,spacing]),np.array([spacing,spacing])))

mem1 = process.memory_info().rss

cplx = ecc.EmbeddedComplex(img*(-1)**dual, 1, input_top_cells=(not dual))
if not dual:
	cplx.impose_upper_star_filtration()
cplx.init_ect(1)

mem2 = process.memory_info().rss

print((mem2 - mem1) / (1024**2))