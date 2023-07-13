import numpy as np

def num_crit_to_spacing(dim, size, n_crit_pts):
	spacings = []
	for ncrit in n_crit_pts:
		spacings.append(int((size/2)*((2/ncrit)**(1/dim))))
	return spacings

def rectangle(sizes,value):
	return np.full(sizes,value)

# def regular_points(sizes, pt_size, pt_spacing, value=1, default_value=0):
# 	arrays = []
# 	for i in range(len(sizes)):
# 		arr = []
# 		for j in range(sizes[i]):
# 			if(j%(pt_spacing[i]+pt_size[i]) < pt_spacing[i]):
# 				arr.append(value)
# 			else:
# 				arr.append(default_value)
# 		arrays.append(np.array([arr]))

# 	res = arrays[0]
# 	for i in range(len(arrays)-1):
# 		res = np.dot(res.T,arrays[i+1])
# 	return res

# def regular_points(sizes, pt_size, pt_spacing, value=1, default_value=0):
# 	arrays = []
# 	for i in range(len(sizes)):
# 		arr = []
# 		for j in range(sizes[i]):
# 			if(j%(pt_spacing[i]+pt_size[i]) < pt_spacing[i]):
# 				arr.append(value)
# 			else:
# 				arr.append(default_value)
# 		arrays.append(np.array(arr))

# 	res = np.full(sizes,1)
# 	for c in np.ndindex(res.shape):
# 		for i in range(len(c)):
# 			res[c] *= arrays[i][c[i]]
# 	return res

def regular_points(sizes, pt_size, pt_spacing, value=1, default_value=0):
	arrays = []
	for i in range(len(sizes)):
		arr = []
		for j in range(sizes[i]):
			if(j%(pt_spacing[i]+pt_size[i]) < pt_spacing[i]):
				arr.append(value)
			else:
				arr.append(default_value)
		arrays.append(np.array(arr))
	a = np.outer(arrays[0],arrays[1])
	for i in range(2, len(sizes)):
		a = np.tensordot(a, arrays[i], axes=0)
	return a
	# if len(sizes)==3:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	return b
	# if len(sizes)==4:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	c = np.tensordot(b,x,axes=0)
	# 	return c
	# if len(sizes)==5:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	c = np.tensordot(b,x,axes=0)
	# 	d = np.tensordot(c,x,axes=0)
	# 	return d
	# if len(sizes)==6:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	c = np.tensordot(b,x,axes=0)
	# 	d = np.tensordot(c,x,axes=0)
	# 	e = np.tensordot(d,x,axes=0)
	# 	return e
	# if len(sizes)==7:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	c = np.tensordot(b,x,axes=0)
	# 	d = np.tensordot(c,x,axes=0)
	# 	e = np.tensordot(d,x,axes=0)
	# 	f = np.tensordot(e,x,axes=0)
	# 	return f
	# if len(sizes)==8:
	# 	x = arrays[0]
	# 	a = np.outer(x,arrays[1])
	# 	b = np.tensordot(a,x,axes=0)
	# 	c = np.tensordot(b,x,axes=0)
	# 	d = np.tensordot(c,x,axes=0)
	# 	e = np.tensordot(d,x,axes=0)
	# 	f = np.tensordot(e,x,axes=0)
	# 	g = np.tensordot(f,x,axes=0)
	# 	return g

def sphere(sizes, radius, thickness, value=1, default_value=0):
	arr = np.full(sizes,default_value)
	center = np.multiply(sizes,0.5)
	for coords in np.ndindex(sizes):
		if abs(np.sum(np.square(np.subtract(coords,center)))**0.5 - radius) < thickness:
			arr[coords] = value
	return arr

def lower_triangle(sizes, value=1, default_value=0):
	arr = np.full(sizes, default_value)
	for coords in np.ndindex(sizes):
		if np.sum(coords) < np.sum(np.multiply(sizes,0.5)):
			arr[coords] = value
	return arr 

def lines(sizes, thickness, spacing, axis=0, value=1, default_value=0):
	arr = np.full(sizes, default_value)
	for coords in np.ndindex(sizes):
		if coords[axis]%(spacing+thickness)+1 > spacing:
			arr[coords] = value
	return arr

