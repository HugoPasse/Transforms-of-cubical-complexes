import numpy as np
import embedded_cubical_complex as ecc


data = np.array([[1,0,0,0],
				[0,1,1,0],
				[1,1,1,1],
				[1,1,0,0],])

cplx = ecc.EmbeddedComplex(data*(-1),input_top_cells=False)

cplx.print_filtration()

print('\n---------------------------------------------')
print('Pre-processing')
cplx.preproc_radon_transform()

directions = [[-1,-1],[-1,1],[1,-1],[1,1]]
for d in directions:
	i = cplx.get_vector_index(d)
	print('\nDirection :',d)
	print('Ordinary critical vertices  :',cplx.get_ordinary_critical_vertices(i))
	print('Ordinary critical values    :',cplx.get_ordinary_critical_values(i))
	print('Classical critical vertices :',cplx.get_classical_critical_vertices(i))
	print('Classical critical values   :',cplx.get_classical_critical_values(i))
print('')

print('\n---------------------------------------------')
direction = [2,2]
print('Computing ect in direction',direction)
dir_id = cplx.get_vector_index(direction)
print('Direction index :',dir_id,end='\n\n')

vertices = cplx.get_classical_critical_vertices(dir_id)
variations = cplx.get_classical_critical_values(dir_id)

print('Classical critical vertices :',vertices)
print('Classical critical values   :',variations)
print('')
dots = []
for i in range(len(vertices)):
	print('Vertex',vertices[i])
	print('  Embedding',end='')
	print(["{0:0.3f}".format(j) for j in cplx.get_vertex_embedding(vertices[i])])
	print('  Scalar product with',direction,':',"{:0.2f}".format(np.dot(cplx.get_vertex_embedding(vertices[i]),direction)))
	print('  Variation',variations[i])
	dots.append(np.around(np.dot(cplx.get_vertex_embedding(vertices[i]),direction),5))
print('')

sort = np.argsort(dots)
print('Sorted dot products :',["{0:0.3f}".format(dots[sort[i]]) for i in sort])
print('Sorted vertices     :',["{:d}".format(vertices[sort[i]]) for i in sort])
print('Sorted variations   :',["{0:0.3f}".format(variations[sort[i]]) for i in sort])

print('\nSumming the variations in sorted order :')
T = [-np.Inf]
Values = [0]
step = 0
for i in sort:
	print("Step",step,":")
	print("  T      =",T)
	print("  Values =",Values)

	if T[-1] == dots[sort[i]]:
		print('  ',dots[sort[i]],'= T[-1]')
		print(' Adding',variations[sort[i]],'to the last element of values')
		Values[-1] += variations[sort[i]]
	else:
		print('  ',dots[sort[i]],'!= T[-1]')
		print('  appening',dots[sort[i]],'to T and',Values[-1] + variations[sort[i]],'to Values')
		T.append(dots[sort[i]])
		Values.append(Values[-1] + variations[sort[i]])
	step += 1
	print('')

print('\n---------------------------------------------')
print('Result :')
ect = cplx.compute_euler_caracteristic_transform(direction)
attributes = ect.get_attributes()
print('T values   :',end='')
print(["{0:0.3f}".format(i) for i in attributes[0]])
print('ECT values :',end='')
print(["{0:0.3f}".format(i) for i in attributes[1]])
