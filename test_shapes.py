import numpy as np

def rectangle(sizes,value):
	return np.full(sizes,value)

def regular_points(sizes, pt_size, pt_spacing, value=1, default_value=0):
	arr = np.full(sizes,default_value)
	for coords in np.ndindex(sizes):
		b = 1
		for i in range(len(coords)):
			if coords[i]%(pt_spacing[i]+pt_size[i]) < pt_spacing[i]:
				b = 0
				break
		if b:
			arr[coords] = 1
	return arr

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

