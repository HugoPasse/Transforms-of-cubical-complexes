import numpy as np

def num_crit_to_spacing(dim, sizes, n_crit_pts):
	spacings = np.zeros((len(sizes),len(n_crit_pts)))
	for i, size in enumerate(sizes):
		for j, ncrit in enumerate(n_crit_pts):
			spacings[i,j] = (size/2)*((2/ncrit)**(1/dim))
	return np.floor(spacings)

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
	return np.outer(arrays[0],arrays[1])

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

