#%%
import numpy as np
import os
import subprocess
from test_shapes import num_crit_to_spacing

expe = 0

# Experiment 0: test
if expe == 0:
	sizes = [10, 20, 30]
	n_dir = 10
	spacings = [1*(i+1) for i in range(3)]
	n_samples = 1
	dual = 1
	dim = 2
	title = 'test'
	transform = 'HT' # 'ECT' or 'Radon' or 'HT'

# Experiment 1: size 
if expe == 1:
	sizes = [20, 40, 100, 500, 1000, 5000]
	n_dir = 100
	spacings = [size//20 for size in sizes]
	n_samples = 10
	dual = 1
	dim = 2
	transform = 'HT' # 'ECT' or 'Radon' or 'HT'
	title = 'our-size-' + transform + '-dim-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

if expe == 2:
	sizes = [1000]
	n_dir = 100
	dim = 2
	n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000]
	spacings = num_crit_to_spacing(dim, sizes[0], n_crit_pts)
	n_samples = 10
	dual = 1
	transform = 'HT' # 'ECT' or 'Radon' or 'HT'
	title = 'our-crit-pts-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

#%% Experiment
overwrite_lock = str(np.random.rand())
path_to_savings = 'timings/' + title + '-' + overwrite_lock
print('######################')
print('transform:', transform)
print('sizes:', sizes)
print('spacings:\n', spacings)
print('nbre directions:', n_dir)
print('dual:', dual)

with open(path_to_savings + '-logs.txt', 'a+') as file:
	file.write('############\n\n')
	file.write(f'transform: {transform} (if HT, kernel=cos)\n')
	file.write(f'sizes: {sizes}\n')
	file.write(f'spacings:\n{spacings}\n')
	file.write(f'nbre directions: {n_dir}\n')
	file.write(f'dual: {dual}\n')
	file.write('(time (s),memory (KB)): \n init cplx \n upper_star \n pre_processing \n transform \n total \n\n')
	file.close()

result = np.zeros((n_samples,len(sizes),len(spacings),5,2))
n_crit_pts = np.zeros((n_samples,len(sizes),len(spacings), 2))
for _ in range(n_samples):
	for i, size in enumerate(sizes):
		for j, spacing in enumerate(spacings):
			print(f'sample = {_} | size = {size} | spacing = {spacing}')
			cmd = 'python3 subprocess_mem.py ' +  str(size) + ' ' + str(spacing) + ' ' + str(n_dir) + ' ' + str(dual) + ' ' + str(dim) + ' ' + path_to_savings + ' ' + transform
			output = subprocess.check_output(cmd, shell=True)
			temp_res = [float(_.decode()) for _ in output.split()]
			result[_,i,j] = np.array([(temp_res[2*_],temp_res[2*_+1]) for _ in range(4)]+[(sum(temp_res[2*_] for _ in range(4)), sum(temp_res[2*_+1] for _ in range(4)))])
			n_crit_pts[_,i,j] = np.array([int(temp_res[-2]), int(temp_res[-1])])
			with open(path_to_savings+'-logs.txt', 'a+') as file:
				file.write(f'\n nbr critical points: {n_crit_pts[_,i,j]}\n')
				file.write(f'result:\n{result[_,i,j]}\n')
				file.close()
np.savez(path_to_savings, result=result, n_crit_pts=n_crit_pts)
print('Results saved.')
# %%
