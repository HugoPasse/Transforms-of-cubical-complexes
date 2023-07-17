#%%
import numpy as np
import os
import subprocess

def main(transform='HT', range_val=256): # transform = 'ECT' or 'Radon' or 'HT'
	size = 2
	n_dir = 100 
	dual = 1
	title = 'our-fashion_mnist-res-' + str(range_val) + '-' + transform + '-n_dir-' + str(n_dir) + '-dual-' + str(dual)

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
	result = np.zeros((5,2))
	n_crit_pts = np.zeros((size, 2))
	total_ext = np.zeros(2)
	for _ in range(size):
		print(f'index = {_}')
		cmd = 'python3 fashion_mnist_our_mem.py ' +  str(size) + ' ' + str(n_dir) + ' ' + str(dual) + ' ' + path_to_savings + ' ' + transform + ' ' + str(range_val) + ' ' + str(_)
		# cmd = '/usr/bin/time --output=' + path_to_savings + '-total-logs.txt -f "%U %M" python3 fashion_mnist_our_mem.py ' +  str(size) + ' ' + str(n_dir) + ' ' + str(dual) + ' ' + path_to_savings + ' ' + transform + ' ' + str(range_val) + ' ' + str(_)
		output = subprocess.check_output(cmd, shell=True)
		temp_res = [float(_.decode()) for _ in output.split()]
		result += np.array([(temp_res[2*_],temp_res[2*_+1]) for _ in range(4)]+[(sum(temp_res[2*_] for _ in range(4)), sum(temp_res[2*_+1] for _ in range(4)))])
		n_crit_pts[_] = np.array([int(temp_res[-2]), int(temp_res[-1])])
		with open(path_to_savings+'-logs.txt', 'a+') as file:
			file.write(f'\n nbr critical points (cla,ord): {n_crit_pts[_]}\n')
			file.write(f'result so far:\n{result}\n')
			file.close()
		# with open(path_to_savings+'-total-logs.txt', 'r') as file:
		# 	line = file.readline().split()
		# 	total_ext += np.array([float(line[0]), int(line[1])])
		# os.remove(path_to_savings + '-total-logs.txt')
	np.savez(path_to_savings, result=result, n_crit_pts=n_crit_pts, total_ext=total_ext)
	print('Results saved in:', path_to_savings)

#%%
main('HT', range_val=2)
main('ECT', range_val=2)
main('Radon', range_val=2)

# %%
