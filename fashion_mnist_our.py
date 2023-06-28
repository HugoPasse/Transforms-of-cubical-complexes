#%%
import time
import numpy as np
clear_line = '\r'+' '*50+'\r'
#%%
import embedded_cubical_complex as ecc
import test_shapes

def timing_our(img, directions, dual=False):
	print('Timing our...', end=clear_line)
	# Construct complex
	tic = time.perf_counter()
	cplx = ecc.EmbeddedComplex(img*(-1)**dual, 1, input_top_cells=(not dual))
	toc = time.perf_counter()
	time_cplx = toc-tic
	
	# Impose upper star
	time_ustar = 0
	if not dual:
		tic = time.perf_counter()
		cplx.impose_upper_star_filtration()
		toc = time.perf_counter()
		time_ustar = toc-tic

	# Pre-processing
	tic = time.perf_counter()
	cplx.init_ect(1)
	toc = time.perf_counter()
	time_preproc = toc-tic

	# Computing ECT
	euler_tr = 0
	tic = time.perf_counter()
	for i in range(len(directions)):
		euler_tr = cplx.compute_euler_caracteristic_transform(directions[i])
	toc = time.perf_counter()
	time_ECT = toc-tic
	
	time_total = time_cplx + time_preproc + time_ECT
	
	print('Done.', end=clear_line)

	return cplx, time_cplx, time_ustar, time_preproc, time_ECT, time_total

#%%
import demeter.euler as euler
import demeter.misc as misc
import demeter.directions as dirs

def timing_demeter(img, directions, T=32):
	print('Timing dem...', end=clear_line)
	# Initializing cubical complex
	tic = time.perf_counter()
	cplx = euler.CubicalComplex(img)
	toc = time.perf_counter()
	time_init = toc-tic

	# "Complexifying"
	tic = time.perf_counter()
	seed = cplx.complexify()
	toc = time.perf_counter()
	time_complex = toc-tic

	# Computing ECT 
	tic = time.perf_counter()
	ect = seed.ECT(directions, T)
	toc = time.perf_counter()
	time_ECT = toc-tic

	time_total = time_init + time_complex + time_ECT
	print('Done.', end=clear_line)
	return cplx, seed, ect, time_init, time_complex, time_ECT, time_total
# %% Fashion_MNIST
#%%
import tensorflow_datasets as tfds

def timing_dataset(dataset, n_dir, T, num_thresholds, time_our=True, time_dem=True, stop='full', dual=False):
	print('########### Timing dataset ###########')
	print('Dataset:', dataset)
	print('Stop:', stop)
	print('Nbre directions:', n_dir)
	print('Num_thresholds:', num_thresholds)

	ds = tfds.as_numpy(tfds.load(dataset, split='train', shuffle_files=True, as_supervised=True))
	directions = np.random.rand(n_dir,2)

	T_our = np.zeros(5)
	crit_pts = []
	
	quantiles = [k/(num_thresholds+1) for k in range(1,num_thresholds+1)]

	overwrite_lock = str(np.random.rand())
	path_to_savings = 'timings/timing-logs-T-{}-num_thresh-{}-ndir-{}-'.format(T,num_thresholds,n_dir)+ dataset + '-' + overwrite_lock
	with open(path_to_savings+'.txt', 'a+') as file:
		file.write('###### Timing ######\n\n')
		file.write('Nbre directions: {}\n'.format(n_dir))
		file.write('T in sec')
		file.write('Our: \t T init cplx \t|\t T pre_processing \t|\t T ECT \t|\t T Total \n\n')
		file.close()

	for i, (img, label) in enumerate(ds):
		img = img.reshape((img.shape[0],img.shape[1]))
		values = np.unique(img)[1:]
		val_quantiles = [np.quantile(values, q) for q in quantiles]

		with open(path_to_savings+'.txt', 'a+') as file:
			file.write('\nIndex = {}\n'.format(i))
			file.write('Label = {}\n'.format(label))
		if stop=='full' or i < int(stop):
			if time_our: 
				img_masked = np.zeros_like(img)
				for val_pix in val_quantiles:
					img_masked[img>val_pix] = val_pix
				our = timing_our(img_masked, directions, dual)
				cplx = our[0]
				T_our += our[1:]
				cplx.init_hybrid_transform(1)
				n_crit = len(cplx.get_critical_vertices(0))
				crit_pts.append(n_crit)
				with open(path_to_savings+'.txt', 'a+') as file:
					file.write('Nbr critical points: {}\n'.format(n_crit))
					file.write('Our: {}\n'.format(our[1:]))
					file.write('Our current total: {}\n'.format(T_our))
					file.close()
		else:
			break

	np.savez(path_to_savings, crit_pts = np.array(crit_pts), timings_our=T_our)
	print('########### Timed ###########')
	return crit_pts, T_our

#%% Test
# n_dir = 10
# dataset = 'fashion_mnist'
# N, T_our = timing_dataset(dataset, n_dir, 100, 10, stop='2')

# %%
# TODO: DONT FORGET TO SCREEN BEFORE LAUNCHING THE TASK

n_dir = 100
dataset = 'fashion_mnist'

# T = 100, num_thresholds = 10
num_thresholds = 10
timing_dataset(dataset, n_dir, 100, num_thresholds, dual=True)