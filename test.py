import embedded_cubical_complex as ecc

import matplotlib.pyplot as plt
import numpy as np
import cv2

from pylab import cm,imshow,colorbar,show

def _2D_array_plot(A,ranges,save_path=None):
    img = imshow(A,extent=ranges,cmap=cm.RdBu)
    colorbar(img)
    if save_path is not None:
        plt.savefig(save_path)
    show()
    

def plot_2D_hybrid_transform(embedded_complex, kernel, start_x, end_x, num_x, start_y, end_y, num_y):
    X = np.linspace(start_x, end_x, num_x)
    Y = np.linspace(start_y, end_y, num_y)
    directions = np.transpose([np.repeat(X, num_y), np.tile(Y, num_x)])
    Z = embedded_complex.compute_hybrid_transform(kernel, directions)
    Z = np.transpose(np.reshape(Z,(num_x,num_y)))
    _2D_array_plot(Z,[start_x,end_x,start_y,end_y])
    return Z

def plot_2D_hybrid_transform_exp_norm(embedded_complex, start_x, end_x, num_x, start_y, end_y, num_y, normalize=False, save_path=None):
    X = np.linspace(start_x, end_x, num_x)
    Y = np.linspace(start_y, end_y, num_y)
    directions = np.transpose([np.repeat(X, num_y), np.tile(Y, num_x)])
    S = embedded_complex.compute_hybrid_transform("sin", directions)
    C = embedded_complex.compute_hybrid_transform("cos", directions)
    Z = np.add(np.square(C),np.square(S))
    if normalize:
        scalar_pdt = np.linalg.norm(directions, axis=1)
        Z = np.divide(Z,scalar_pdt)
    Z = np.flip(np.transpose(np.reshape(Z, (num_x,num_y))), axis=0)
    _2D_array_plot(Z,[start_x,end_x,start_y,end_y],save_path)
    return Z

def embedded_cubical_complex_from_path(path):
    img = cv2.imread(path)
    img = np.multiply(img,-1)
    cplx = ecc.EmbeddedComplex(img)
    return cplx

def plot_radon_transform(cplx, direction, title=None):
    radon = cplx.compute_radon_transform(direction)
    attributes = radon.get_attributes()
    print(attributes)
    if(len(attributes[0]) > 0):
        plt.plot([-1,attributes[0][0]],[0,0],'b')
        plt.plot([attributes[0][len(attributes[0])-1],1],[0,0],'b')
        for i in range(len(attributes[0])-1):
            plt.plot([attributes[0][i], attributes[0][i+1]], [attributes[1][i], attributes[1][i]], 'b')
        for i in range(len(attributes[2])):
            plt.plot([attributes[2][i]], [attributes[3][i]], 'b', marker='o')
    else:
        plt.plot([-1,1],[0,0])
    if(title is not None):
        plt.title(title)
    plt.show()

def plot_euler_caracteristic_transform(embedded_complex, direction, title=None):
    ect = cplx.compute_euler_caracteristic_transform(direction)
    attributes = ect.get_attributes()
    if(len(attributes[0]) > 0):
        plt.plot([-1,attributes[0][0]],[0,0],'b')
        plt.plot([attributes[0][len(attributes[0])-1],1], [attributes[1][len(attributes[1])-1],attributes[1][len(attributes[1])-1]], 'b')
        for i in range(len(attributes[0])-1):
            plt.plot([attributes[0][i], attributes[0][i+1]], [attributes[1][i], attributes[1][i]], 'b')
    else:
        plt.plot([-1,1],[0,0])
    if(title is not None):
        plt.title(title)
    plt.show()

"""
data = np.array([[-1,0,0],[0,-1,0],[0,0,-1]])
direction = [1,1]
cplx = ecc.EmbeddedComplex(data)
cplx.print_filtration()
plot_2D_hybrid_transform_exp_norm(cplx, -5,5,11,-5,5,11)

plt.waitforbuttonpress()
"""

direction = [1,1]
np.set_printoptions(threshold=np.inf)
#img = cv2.imread("/home/hugo/Documents/L3/S2/stage/implementation/images/belgian-fieldstone.png",cv2.IMREAD_GRAYSCALE)
img = np.array([[0,0,0,0],[0,0,1,0],[0,1,0,1],[0,0,1,0]])
#img = np.multiply(img,1/255.0)
#img = np.flip(np.flip(img), 1)

cplx = ecc.EmbeddedComplex(img,1)
cplx.print_filtration()
#cplx.print_embedding()
#cplx.print_critical_vertices()
#cplx.print_critical_multiplicity()
cplx.init_radon_transform()
plot_radon_transform(cplx,direction,"Radon")

#plot_euler_caracteristic_transform(cplx, direction);
#plot_2D_hybrid_transform_exp_norm(cplx,-100,100,1000,-100,100,1000, normalize=True)

#plot_radon_transform(cplx,direction)

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