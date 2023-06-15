#%%
import time
import numpy as np

#%%
import embedded_cubical_complex as ecc
import test_shapes

def timing_our(img, directions):
	print('Timing our...', end='\r')
	# Construct complex
	tic = time.perf_counter()
	cplx = ecc.EmbeddedComplex(img,1)
	toc = time.perf_counter()
	time_cplx = toc-tic

	# Pre-processing
	tic = time.perf_counter()
	cplx.init_hybrid_transform(1)
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
	
	print('Done.', end='\r')

	return cplx, time_cplx, time_preproc, time_ECT, time_total

#%%
import demeter.euler as euler
import demeter.misc as misc
import demeter.directions as dirs

def timing_demeter(img, directions, T=32):
	print('Timing dem...', end='\r')
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
	print('Done.', end='\r')
	return cplx, seed, ect, time_init, time_complex, time_ECT, time_total

#%%
def num_crit_to_spacing(size, dim, ncrit):
	p = (size/2)*((2/float(ncrit))**(1/float(dim)))
	return int(p)

#%%
def timings_by_critical_points(n, n_crit_pts, directions, T=32, title = '', n_samples=1):
	print('Timing by critical points...')
	T_our = np.zeros((n_samples,len(n_crit_pts),4))
	T_dem = np.zeros((n_samples,len(n_crit_pts),4))
	N = np.zeros(len(n_crit_pts))
	with open('timings/timing_by_critical_points-'+ title + '-' + str(np.random.rand())+'.txt', 'a+') as file:
		file.write('Timing by critical points\n\n')
		file.write('Our time: \t init cplx | pre_processing | ECT | Total\n')
		file.write('Demeter time (T={}): \t init | complexifying | ECT | Total\n'.format(T))
		for j in range(len(n_crit_pts)):
			print('Nbr pts critiques: %s (computing...)' %n_crit_pts[j])
			for i in range(n_samples):
				p = num_crit_to_spacing(n, 2, n_crit_pts[j])
				img = test_shapes.regular_points((n,n),np.array([p,p]),np.array([p,p]))
				our = timing_our(img, directions)
				cplx, T_our[i,j,:] = our[0], our[1:]
				N[j] = len(cplx.get_critical_vertices(0))
				file.write('Nbr critical points: {}\n'.format(N[j]))
				file.write('Our: {}\n'.format(our[1:]))
				T_dem[i,j,:] = timing_demeter(img, directions, T = T)[3:]		
				file.write('Dem: {}\n\n'.format(T_dem[i,j,:]))
		file.close()
	print('Timed by critical points.')
	return N, T_our, T_dem

#%%
def timings_by_size(sizes, spacings, directions, T=32, title='', n_samples=1):
	print('Timing by sizes...')
	T_our = np.zeros((n_samples, len(sizes), 4))
	T_dem = np.zeros((n_samples, len(sizes), 4))
	N = np.zeros(len(sizes))
	with open('timings/timing_by_size-'+ title + '-' + str(np.random.rand())+'.txt', 'a+') as file:
		file.write('Timing by size\n\n')
		file.write('Sizes:' + str(sizes) + '\n\n')
		file.write('Our time: init cplx | pre_processing | ECT | Total\n')
		file.write('Demeter time (T={}): init | complexifying | ECT | Total\n'.format(T))
		for j in range(len(sizes)):
			print('Taille de l\'image: %s (computing...)' %sizes[j])
			p = spacings[j]
			n = sizes[j]
			for i in range(n_samples):
				img = test_shapes.regular_points((n,n), np.array([p,p]),np.array([p,p]))
				our = timing_our(img, directions)
				cplx, T_our[i,j,:] = our[0], our[1:]
				N[j] = len(cplx.get_critical_vertices(0))
				file.write('Taille de l\'image: {}\n'.format(sizes[j]))
				file.write('Nbr critical points: {}\n'.format(N[j]))
				file.write('Our: {}\n'.format(our[1:]))
				T_dem[i,j,:] = timing_demeter(img, directions, T = T)[3:]
				file.write('Dem: {}\n\n'.format(T_dem[i,j,:]))
		file.close()
	print('Timed by sizes.')
	return N, T_our, T_dem

# %%
directions = np.random.rand(100,2)

#%%
N, T_our, T_dem = timings_by_critical_points(100, [10, 20, 30], directions)

# %%
sizes = [10**k for k in range(2, 4)]
spacings = [int(size/20) for size in sizes]
N, T_our, T_dem = timings_by_size(sizes, spacings, directions, T=100)
# %%
