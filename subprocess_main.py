#%%
import numpy as np
import os
import subprocess
from test_shapes import num_crit_to_spacing

def main(expe=0, transform='HT'): # transform = 'ECT' or 'Radon' or 'HT'
	# Experiment 0: test
	if expe == 0:
		sizes = [40, 60, 80]
		n_dir = 10
		spacings = [size//20 for size in sizes]
		n_samples = 1
		dual = 1
		dim = 2
		title = 'test'

	# Experiment 1: size 
	if expe == 1:
		sizes = [40, 100, 500, 1000, 5000]
		n_dir = 100
		spacings = [size//20 for size in sizes]
		n_samples = 10
		dual = 1
		dim = 2
		title = 'our-size-' + transform + '-n-dir-' + str(n_dir) + '-n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment 2: critical points
	if expe == 2:
		sizes = [1000]
		n_dir = 100
		dim = 2
		n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000]
		spacings = num_crit_to_spacing(dim, sizes[0], n_crit_pts)
		n_samples = 10
		dual = 1
		title = 'our-crit-pts-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment 3,4: dimension
	if expe==3:
		sizes = [292]
		n_dir = 100
		dim = 3
		spacings = [14]
		n_samples = 10
		dual = 1
		title = 'our-dim-3-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)
	
	if expe==4:
		sizes = [70]
		n_dir = 100 
		dim = 4
		spacings = [3]
		n_samples = 10
		dual = 1
		title = 'our-dim-4-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	if expe==5:
		sizes = [30]
		n_dir = 100 
		dim = 5
		spacings = [1]
		n_samples = 10
		dual = 1
		title = 'our-dim-5-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)
	
	if expe==6:
		sizes = [17]
		n_dir = 100 
		dim = 6
		spacings = [1]
		n_samples = 10
		dual = 1
		title = 'our-dim-6-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment
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
		file.write('(time (s),memory (KB)): \n init cplx \n upper_star \n pre_processing \n transform \n total \n ext total \n\n')
		file.close()
	m = max(len(sizes),len(spacings))
	result = np.zeros((n_samples,m,5,2))
	n_crit_pts = np.zeros((n_samples,m, 2))
	total_ext = np.zeros((n_samples,m, 2))
	for _ in range(n_samples):
		for i in range(m):
			size = sizes[min(len(sizes)-1,i)]
			spacing = spacings[min(len(spacings)-1,i)]
			print(f'sample = {_} | size = {size} | spacing = {spacing}')
			cmd = '/usr/bin/time --output=' + path_to_savings + '-total-logs.txt -f "%U %M" python3 subprocess_mem.py ' +  str(size) + ' ' + str(spacing) + ' ' + str(n_dir) + ' ' + str(dual) + ' ' + str(dim) + ' ' + path_to_savings + ' ' + transform
			output = subprocess.check_output(cmd, shell=True)
			temp_res = [float(_.decode()) for _ in output.split()]
			result[_,i] = np.array([(temp_res[2*_],temp_res[2*_+1]) for _ in range(4)]+[(sum(temp_res[2*_] for _ in range(4)), sum(temp_res[2*_+1] for _ in range(4)))])
			n_crit_pts[_,i] = np.array([int(temp_res[-2]), int(temp_res[-1])])
			with open(path_to_savings+'-logs.txt', 'a+') as file:
				file.write(f'\n nbr critical points (cla,ord): {n_crit_pts[_,i]}\n')
				file.write(f'result:\n{result[_,i]}\n')
				file.close()
			with open(path_to_savings+'-total-logs.txt', 'r') as file:
				line = file.readline().split()
				total_ext[_,i] = np.array([float(line[0]), int(line[1])])
	np.savez(path_to_savings, result=result, n_crit_pts=n_crit_pts, total_ext=total_ext)
	print('Results saved in:', path_to_savings)

#%%
# for expe in [0,1,2]:
# 	for transform in ['HT', 'Radon', 'ECT']:
# 		main(expe,transform)

#%%
for transform in ['HT', 'Radon', 'ECT']:
	main(3,transform)
	main(4,transform)
	main(5,transform)
	main(6,transform)