import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
import time
import sys

sys.path.append("/home/hugo/Documents/L3/S2/stage/demeter")

import embedded_cubical_complex as ecc

import demeter.euler as euler
import demeter.misc as misc
import demeter.directions as dirs

import test_shapes

i = 10
dim = 2
img = test_shapes.regular_points((2000,2000),np.array([i,i]),np.array([i,i]))

# file = 'test_images/seed.tif'
# img = tf.imread(file)
# img[img > 0] = 1
# print(img.shape)

ndirs = 10
directions = np.random.rand(ndirs,dim)

##################################################

# print("---- Demeter ----")

# print("-- Construct complex --")
# tic = time.perf_counter()
# cplx = euler.CubicalComplex(img).complexify()
# max_vox = misc.find_tip(cplx.cells[0], 2,1,0)
# # coords,_, _, _,_ = misc.rotateSVD(cplx.cells[0], max_vox)
# toc = time.perf_counter()
# print(toc-tic)

# print("-- Compute " + str(len(directions)) + " directions --")
# tic = time.perf_counter()
# ect = cplx.ECT(directions, T=32, verts=cplx.cells[0], bbox=(-25,25))
# toc = time.perf_counter()
# print(toc-tic)

##################################################

print("\n---- Our implementation ----")

print("-- Construct complex --")
tic = time.perf_counter()
cplx = ecc.EmbeddedComplex(img,1)
toc = time.perf_counter()
print(toc-tic)

print("-- Init --")
tic = time.perf_counter()
cplx.init_ect()	# DANGER HERE IN PARALLEL
toc = time.perf_counter()
print(toc-tic)

print("-- Compute " + str(len(directions)) + " directions --")
tic = time.perf_counter()
for d in directions:
	cplx.compute_euler_caracteristic_transform(d)
toc = time.perf_counter()
print(toc-tic)

##################################################