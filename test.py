import embedded_cubical_complex as ecc

import matplotlib.pyplot as plt
import numpy as np
import cv2

from pylab import cm,imshow,colorbar,show

def _2D_array_plot(A,ranges):
    img = imshow(A,extent=ranges,cmap=cm.RdBu)
    colorbar(img)
    show()

def plot_2D_hybrid_transform(embedded_complex, kernel, start_x, end_x, num_x, start_y, end_y, num_y):
    X = np.linspace(start_x, end_x, num_x)
    Y = np.linspace(start_y, end_y, num_y)
    directions = np.transpose([np.repeat(X, num_y), np.tile(Y, num_x)])
    Z = embedded_complex.compute_hybrid_transform(kernel, directions)
    Z = np.transpose(np.reshape(Z,(num_x,num_y)))
    _2D_array_plot(Z,[start_x,end_x,start_y,end_y])
    return Z

def plot_2D_hybrid_transform_exp_norm(embedded_complex, start_x, end_x, num_x, start_y, end_y, num_y):
    X = np.linspace(start_x, end_x, num_x)
    Y = np.linspace(start_y, end_y, num_y)
    directions = np.transpose([np.repeat(X, num_y), np.tile(Y, num_x)])
    S = embedded_complex.compute_hybrid_transform("sin", directions)
    C = embedded_complex.compute_hybrid_transform("cos", directions)
    Z = np.transpose(np.reshape(np.add(np.square(C),np.square(S)),(num,num)))
    _2D_array_plot(Z,[start_x,end_x,start_y,end_y])
    return Z

def embedded_cubical_complex_from_path(path):
    img = cv2.imread(path)
    img = np.multiply(img,-1)
    cplx = ecc.EmbeddedComplex(img)
    return cplx


"""
img = cv2.imread("/home/hugo/Documents/L3/S2/stage/implementation/images/belgian-fieldstone.png")
img = np.multiply(img,-1)
cplx = ecc.EmbeddedComplex(img)

start = -50
end = 50
num = 100

x = np.linspace(start, end, num)
directions = np.transpose([np.repeat(x, num), np.tile(x, num)])

C = cplx.compute_hybrid_transform("cos",directions)
S = cplx.compute_hybrid_transform("sin",directions)
Z = np.transpose(np.reshape(np.add(np.square(C),np.square(S)),(num,num)))
_2D_array_plot(Z,[start,end,start,end])
"""