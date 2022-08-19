import numpy as np
from matplotlib import pyplot
from matplotlib import image
import os
import time

import embedded_cubical_complex as embedded_cc


def load_data(path):
	images = []
	for subdir in os.listdir(path):
		subdir_path = os.path.join(path, subdir)
		for filename in os.listdir(subdir_path):
			filename_path = os.path.join(subdir_path, filename)
			if os.path.isfile(filename_path):
				img = image.imread(filename_path)
				img = np.multiply(img,255).astype(int)
				images.append(img)
	return images

def compute_transform_on_img(img, vectors):
	cplx = embedded_cc.EmbeddedComplex(img)
	cplx.compute_hybrid_transform("exp", vectors)

def random_vectors(number, dimension, val_range):
	return (np.multiply(np.random.rand(number,dimension), val_range))

def compute_on_dataset(path,num=100,val_range=10):
	images = load_data(path)
	vectors = random_vectors(num,2,val_range)
	print("Image count : " + str(len(images)))
	start = time.time()
	for img in images:
		compute_transform_on_img(img,vectors)
	end = time.time()
	dt = end - start
	print("Time : " + str(dt) + " s")


compute_on_dataset('/home/hugo/Documents/L3/S2/stage/Transforms-of-cubical-complexes/tests/python/mnist')