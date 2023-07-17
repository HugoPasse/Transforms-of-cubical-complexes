import numpy as np
import time
import matplotlib.pyplot as plt

import embedded_cubical_complex as ecc
from test_shapes import regular_points

ncores = 8 
ncplx = 100
sizes = [50,100,250,500,2000]


cores = [i+1 for i in range(ncores)]

TS = []
for size in sizes:
	T = []
	for nc in cores:
		print('size : ',size,'cores :',nc)
		t = 0
		for i in range(ncplx): 
			data = np.array(regular_points(tuple([size]*2),np.array([1]*2),np.array([1]*2)))
			cplx = ecc.EmbeddedComplex(data, 1, input_top_cells=True)
			tic = time.perf_counter()
			cplx.preproc_radon_transform(num_jobs=nc)
			toc = time.perf_counter()
			t += toc-tic
		T.append(t/ncplx)
	TS.append(T)

for i in range(len(sizes)):
	plt.plot(cores,np.array(TS[i])/TS[i][0],label='Size = '+str(sizes[i]))

plt.legend()
plt.xlabel('Number of cores used')
plt.ylabel('Time / time with one core')
plt.title('Parallelization performance on random binary complexes')

plt.savefig("parallel_pref")

plt.show()