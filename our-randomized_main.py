#%%
import numpy as np
import os
import subprocess

def main(range_val, transform='HT'): # transform = 'ECT' or 'Radon' or 'HT'
	size = 10 # What should the size be?
	# size = 1000 # What should the size be? 
	n_dir = 100 
	dim = 2
	n_samples = 10
	dual = 1
	title = 'our-range-values-' + str(range_val) + '-' + transform + '-n_dir-' + str(n_dir) + 'n-samples-' + str(n_samples) + '-dual-' + str(dual) + '-dim-' + str(dim)

	# Experiment
	overwrite_lock = str(np.random.rand())
	path_to_savings = 'timings/' + title + '-' + overwrite_lock
	print('######################')
	print('transform:', transform)
	print('size:', size)
	print('nbre directions:', n_dir)
	print('dual:', dual)
	print('range of values of img:', range_val)

	with open(path_to_savings + '-logs.txt', 'a+') as file:
		file.write('############\n\n')
		file.write(f'transform: {transform} (if HT, kernel=cos)\n')
		file.write(f'range_val: {range_val}\n')
		file.write(f'size: {size}\n')
		file.write(f'nbre directions: {n_dir}\n')
		file.write(f'dual: {dual}\n')
		file.write('(time (s),memory (KB)): \n init cplx \n upper_star \n pre_processing \n transform \n total\n\n')
		file.close()
	result = np.zeros((n_samples,5,2))
	n_crit_pts = np.zeros((n_samples,2))
	total_ext = np.zeros((n_samples,2))
	for _ in range(n_samples):
		print(f'sample = {_} | size = {size}')
		cmd = '/usr/bin/time --output=' + path_to_savings + '-total-logs.txt -f "%U %M" python3 our-randomized_mem.py ' +  str(size) + ' ' + str(n_dir) + ' ' + str(dual) + ' ' + str(dim) + ' ' + path_to_savings + ' ' + transform + ' ' + str(range_val)
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
for transform in ['HT', 'Radon', 'ECT']:
	for range_val in [2, 5, 10]:
		main(range_val, transform)
# for transform in ['HT', 'Radon', 'ECT']:
# 	for range_val in [2, 5, 10, 50, 100, 200, 300, 500, 1000]:
# 		main(range_val, transform)

# TODO: test beforehand