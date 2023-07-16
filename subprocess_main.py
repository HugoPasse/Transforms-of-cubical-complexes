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
		sizes = [2000]
		n_dir = 100
		dim = 2
		spacings = [400, 300, 250, 200, 150, 125, 110, 100, 90, 80, 70, 60, 50, 45, 42, 40, 35, 30, 25, 20, 15, 10, 9, 8, 7, 6, 5, 4, 3, 2] # n_crit_pts = [4, 9, 16, 25, 36, 64, 81, 100, 121, 144, 196, 256, 400, 484, 529, 625, 784, 1089, 1600, 2500, 4356, 10000, 12321, 15625, 20164, 27556, 40000, 62500, 110889, 250000]
		n_samples = 10
		dual = 1
		title = 'our-crit-pts-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# # Experiment 2: critical points (old)
	# if expe == 2:
	# 	sizes = [1000]
	# 	n_dir = 100
	# 	dim = 2
	# 	n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000]
	# 	spacings = num_crit_to_spacing(dim, sizes[0], n_crit_pts)
	# 	n_samples = 10
	# 	dual = 1
	# 	title = 'our-crit-pts-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

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
	
	# Experiment 7: size dim 3
	if expe == 7:
		sizes = [40, 100, 200, 400, 500, 1000, 5000]
		n_dir = 100
		spacings = [size//20 for size in sizes]
		n_samples = 10
		dual = 1
		dim = 3
		title = 'our-size-' + transform + '-n-dir-' + str(n_dir) + '-n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment 8: critical points dim 3
	if expe == 8:
		sizes = [500]
		n_dir = 100
		dim = 3
		spacings = [1, 2, 3, 4, 5, 6, 7, 10, 20, 30, 50, 60, 70, 100, 200] # n_crit_pts = [15625000, 1953125, 571787, 238328, 125000, 68921, 42875, 15625, 1728, 512, 125, 64, 27, 8, 1]
		n_samples = 10
		dual = 1
		title = 'our-crit-pts-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment
	overwrite_lock = str(np.random.rand())
	path_to_savings = 'timings/' + title + '-' + overwrite_lock
	print('######################')
	print('transform:', transform)
	print('sizes:', sizes)
	print('spacings:\n', spacings)
	print('nbre directions:', n_dir)
	print('dual:', dual)
	print('dim:', dim)

	with open(path_to_savings + '-logs.txt', 'a+') as file:
		file.write('############\n\n')
		file.write(f'transform: {transform} (if HT, kernel=cos)\n')
		file.write(f'sizes: {sizes}\n')
		file.write(f'spacings:\n{spacings}\n')
		file.write(f'nbre directions: {n_dir}\n')
		file.write(f'dual: {dual}\n')
		file.write(f'dim: {dim}\n')
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
	os.remove(path_to_savings + '-total-logs.txt')
	np.savez(path_to_savings, result=result, n_crit_pts=n_crit_pts, total_ext=total_ext)
	print('Results saved in:', path_to_savings)

#%%
# for expe in [0,1,2]:
# 	for transform in ['HT', 'Radon', 'ECT']:
# 		main(expe,transform)

#%%
# for transform in ['HT', 'Radon', 'ECT']:
# 	main(3,transform)
# 	main(4,transform)
# 	main(5,transform)
# 	main(6,transform)

#%%
for transform in ['HT', 'Radon', 'ECT']:
	main(2,transform)