#%%
import time
import numpy as np
clear_line = '\r'+' '*50+'\r'
#%%
import embedded_cubical_complex as ecc
import test_shapes

def timing_our(img, directions):
	print('Timing our...', end=clear_line)
	# Construct complex
	tic = time.perf_counter()
	cplx = ecc.EmbeddedComplex(img,1)
	toc = time.perf_counter()
	time_cplx = toc-tic

	# Pre-processing
	tic = time.perf_counter()
	cplx.init_ect(1)
	# cplx.init_hybrid_transform(1)
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

	return cplx, time_cplx, time_preproc, time_ECT, time_total

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

#%%
def num_crit_to_spacing(dim, sizes, n_crit_pts):
	spacings = np.zeros((len(sizes),len(n_crit_pts)))
	for i, size in enumerate(sizes):
		for j, ncrit in enumerate(n_crit_pts):
			spacings[i,j] = (size/2)*((2/ncrit)**(1/dim))
	return np.floor(spacings)

#%%
def timing(sizes, spacings, directions, Ts, title='', n_samples=1, time_our=True, time_dem=True, negative=False):
	print('########### Timing ###########')
	print('sizes:', sizes)
	print('spacings:\n', spacings)
	print('Nbre directions:', len(directions))
	print('Ts (for demeter):', Ts, end='\n'+'-'*50+'\n')

	overwrite_lock = str(np.random.rand())
	path_to_savings = 'timings/timing-'+ title + '-' + overwrite_lock
	with open(path_to_savings+'.txt', 'a+') as file:
		file.write('###### Timing ######\n\n')
		file.write('sizes: {}\n'.format(sizes))
		file.write('spacings:\n{}\n'.format(spacings))
		file.write('Nbre directions: {}\n'.format(len(directions)))
		file.write('Ts (for demeter): {}\n\n'.format(Ts))
		file.write('Our time: \t init cplx \t|\t pre_processing \t|\t ECT \t|\t Total\n')
		file.write('Demeter time: \t init \t\t|\t complexifying \t\t|\t ECT \t|\t Total\n\n')
		file.close()

	T_our = np.zeros((n_samples,len(sizes),len(spacings[0]),4))
	T_dem = np.zeros((n_samples,len(sizes),len(spacings[0]),len(Ts), 4))
	N = np.zeros(len(n_crit_pts))
	for i, size in enumerate(sizes):
		for j in range(len(spacings[0])):
			spacing = spacings[0,j]
			for _ in range(n_samples):
				print('Sample: {}'.format(_+1), end=clear_line)
				raw_img = test_shapes.regular_points((size,size),np.array([spacing,spacing]),np.array([spacing,spacing]))
				img = raw_img if not negative else 1-raw_img
				with open(path_to_savings+'.txt', 'a+') as file:
					file.write('\nSize: {}\n'.format(size))
					file.write('Spacing: {}\n'.format(spacing))
					file.write('Sample: {}\n'.format(_+1))
					file.close()
				if time_our:
					our = timing_our(img, directions)
					cplx, T_our[_,i,j,:] = our[0], our[1:]
					cplx.init_hybrid_transform(1)
					N[j] = len(cplx.get_critical_vertices(0))
					with open(path_to_savings+'.txt', 'a+') as file:
						file.write('Nbr critical points: {}\n'.format(N[j]))
						file.write('Our: {}\n'.format(our[1:]))
						file.close()
				if time_dem:
					for k, T in enumerate(Ts):
						T_dem[_,i,j,k,:] = timing_demeter(img, directions, T = T)[3:]		
						with open(path_to_savings+'.txt', 'a+') as file:
							file.write('Dem (T={}): {}\n'.format(T, T_dem[_,i,j,k,:]))
							file.close()
	np.savez(path_to_savings, timings_our=T_our, timings_dem=T_dem)
	print('########### Timed ###########')
	return N, T_our, T_dem
#%% Test
# dim = 2
# sizes = [20]
# n_crit_pts = [10, 50]
# spacings = num_crit_to_spacing(dim, sizes, n_crit_pts)
# # spacings = np.array([int(size/20) for size in sizes]).reshape((1,len(sizes)))
# n_dir = 50
# directions = np.random.rand(n_dir, dim)
# Ts = [5, 10]
# title = 'test'
# n_samples = 2

# N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=True, time_dem=True)

#%% Timing

# TODO: DONT FORGET TO SCREEN BEFORE LAUNCHING THE TASK

dim = 2
n_samples = 10
n_dirs = [50, 100, 500, 1000]

# Demeter with respect to critical points and T, size = 100
for n_dir in n_dirs:
	title = 'critical_pts_demeter-ndir-'+str(n_dir)+'size-100'
	sizes = [100]
	n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000]
	spacings = num_crit_to_spacing(dim, sizes, n_crit_pts)
	directions = np.random.rand(n_dir,2)
	Ts = [50, 100, 500, 1000]
	N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=False, time_dem=True)

# Demeter with respect to critical points and T, size = 1000
for n_dir in n_dirs:
	title = 'critical_pts_demeter-ndir-'+str(n_dir)+'size-1000'
	sizes = [1000]
	n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000]
	spacings = num_crit_to_spacing(dim, sizes, n_crit_pts)
	Ts = [50, 100, 500, 1000]
	directions = np.random.rand(n_dir,2)
	N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=False, time_dem=True)

# Demeter with respect to sizes and T, crit = 200
for n_dir in [50, 100]:
	title = 'size_demeter_ndir-'+str(n_dir)+'ncrit-200'
	sizes = [20, 40, 100, 500, 1000, 5000]
	spacings = np.array([int(size/20) for size in sizes]).reshape((1,len(sizes))) # 200 critical points
	Ts = [50, 100, 500, 1000]
	directions = np.random.rand(n_dir,2)
	N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=False, time_dem=True)

# Both with respect to critical points and T, n_dir = 1000, size = 100
n_dir = 1000
title = 'critical_pts_ndir-'+str(n_dir)+'size-100'
sizes = [100]
n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000]
spacings = num_crit_to_spacing(dim, sizes, n_crit_pts)
Ts = [32, 50, 100, 500, 1000]
directions = np.random.rand(n_dir,2)
N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=True, time_dem=True)

# Both with respect to critical points and T, n_dir = 1000, size = 1000
title = 'critical_pts_ndir-'+str(n_dir)+'size-1000'
sizes = [1000]
n_crit_pts = [10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000]
spacings = num_crit_to_spacing(dim, sizes, n_crit_pts)
Ts = [32, 50, 100, 500, 1000]
directions = np.random.rand(n_dir,2)
N, T_our, T_dem = timing(sizes, spacings, directions, Ts, title, n_samples, time_our=True, time_dem=True)