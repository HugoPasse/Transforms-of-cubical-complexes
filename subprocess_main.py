import numpy as np
import os
import subprocess

spacings = [1*(i+1) for i in range(10)]
n = 1000
dual = 0

for s in spacings:
	print('n = ',n,'| spacing =',s)
	cmd = 'python3 subprocess_mem.py ' +  str(n) + ' ' + str(s) + ' ' + str(dual)
	os.system(cmd)