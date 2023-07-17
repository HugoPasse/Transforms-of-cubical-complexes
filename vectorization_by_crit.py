import numpy as np
import time
import matplotlib.pyplot as plt

import embedded_cubical_complex as ecc
from test_shapes import regular_points

dim = 2
n = 2500
nv = 10

N_dir = 100

spacings = np.arange(5,25)*10

Tpts = []
Cpts = []

for s in spacings:
	timer = 0
	data = regular_points([n,n],[s,s],[1,1])
	cplx = ecc.EmbeddedComplex(data,input_top_cells=True)
	cplx.preproc_ect()

	for i in range(N_dir):
		direction = np.abs(np.random.rand(dim))
		ect = cplx.compute_euler_caracteristic_transform(direction)
		tic = time.time()
		ect.vectorize(-2,2,nv)
		toc = time.time()
		timer += toc - tic

	Tpts.append(timer / N_dir)
	Cpts.append(len(cplx.get_classical_critical_vertices(3)))

print(Cpts)
print(Tpts)
# plt.plot(Cpts,Tpts)
# plt.xlabel("Number of critical points")
# plt.ylabel("Time")

# # title = "vectorization/vectorization_by_crit_n_" + str(n) + "_dim_" + str(dim) + "_ndir_" + str(N_dir) + "_vectSize_" + str(nv) + ".png"
# # plt.savefig(title)
# plt.show()