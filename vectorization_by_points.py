import numpy as np
import time
import matplotlib.pyplot as plt

import embedded_cubical_complex as ecc
from test_shapes import regular_points

dim = 2
n = 100

data = regular_points([n,n],[1,1],[1,1])
cplx = ecc.EmbeddedComplex(data,input_top_cells=True)
cplx.preproc_ect()
ncrit = len(cplx.get_classical_critical_vertices(3))

N_dir = 100
numpts = np.arange(1,100)*100
Tpts = []

for nv in numpts:
	timer = 0
	print('Vector size :',nv)
	for i in range(N_dir):
		direction = np.abs(np.random.rand(dim))
		ect = cplx.compute_euler_caracteristic_transform(direction)
		tic = time.time()
		ect.vectorize(-2,2,nv)
		toc = time.time()
		timer += toc - tic
	Tpts.append(timer / N_dir)


print(numpts)
print(Tpts)

# plt.plot(numpts,Tpts)

# plt.xlabel('Vector Size')
# plt.ylabel('Time')

# title = "vectorization/vectorization_by_points_n_" + str(n) + "_dim_" + str(dim) + "_ndir_" + str(N_dir) + "_ncrit_" + str(ncrit) +".png"
# print(title)
# plt.savefig(title)

# plt.show()
