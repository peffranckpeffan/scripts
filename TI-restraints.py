import sys
import subprocess
import os
import time
from shutil import copyfile
import numpy as np
import math
import os, errno, glob, shutil, json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#csv xls ods
def intSimps(a, b, n, step):
	h = (b - a) / n
	s = 0
	i = 0
	s = f(a) + f(b)

	for i in np.arange(1,n,step):
		s+= 4*f(a + i*h)
	for i in np.arange(2,n,step):
		s+= 2*f(a + i*h)
	
	s = (s*h)/3
	return s

dir_analyse = '../analysis/'
fold_analyse = sys.argv[1]

num_win = 0
for folder in os.listdir(dir_analyse):
	if os.path.isdir(dir_analyse+folder) and fold_analyse in folder and folder.split(fold_analyse)[1].isdigit():
				num_win = num_win + 1

array_len = len(np.genfromtxt(dir_analyse+fold_analyse+'00/restraints.dat', usecols=(2)).transpose())

simul_data = np.zeros([num_win, 2, array_len], np.float64)
potentials = np.zeros([num_win, array_len], np.float64)
constants = np.zeros([2, num_win], np.float64)
centers = np.zeros([2, num_win], np.float64)
mean_pot = np.zeros([1, num_win], np.float64)


if fold_analyse == 'b':
	rest_info = [['rmsd1'], [2]]
elif fold_analyse == 'l':
	rest_info = [['rmsd2'], [1]]
elif fold_analyse == 'o':
	rest_info = [['orient2'], [2]]
elif fold_analyse == 'p':
	rest_info = [['rmsd1'], [1]]
elif fold_analyse == 'r':
	rest_info = [['rmsd2'], [1]]
elif fold_analyse == 't':
	rest_info = [['posit2', 'posit3'], [1,2]]
elif fold_analyse == 'u':
	rest_info = [['posit3'], [2]]
else:
	sys.exit('Not sure about restraint type!')


line_ref = 0
for window in range(num_win):
	for i in range(len(rest_info[1])):
		simul_data[window, i, 0:array_len] = np.genfromtxt(dir_analyse+fold_analyse+str(window).zfill(2)+'/restraints.dat', usecols=(rest_info[1][i])).transpose()

		with open(dir_analyse+fold_analyse+str(window).zfill(2)+'/colvar.in') as file:
			lines = file.readlines()
			for line in lines:
				if rest_info[0][i] in line:
					line_ref = lines.index(line)

			if fold_analyse != 'o':
				centers[i, window] = float(lines[line_ref+2].split()[1])
			else:
				centers[i, window] = float(lines[line_ref+2].split()[1][1:-1])
			constants[i, window] = float(lines[line_ref+3].split()[1])
	
	if fold_analyse in ['b', 'l', 'p', 'r']:
		potentials[window, 0:array_len] = (simul_data[window, 0, 0:array_len]-centers[0,window])**2
	elif fold_analyse == 'o':
		potentials[window, 0:array_len] = (np.arccos(simul_data[window, i, 0:array_len]))**2
	elif fold_analyse == 't':
		potentials[window, 0:array_len] = (simul_data[window, 0, 0:array_len]-centers[0,window])**2 + (simul_data[window, 1, 0:array_len]-centers[1,window])**2
	elif fold_analyse == 'u':
		potentials[window, 0:array_len] = constants[0, window]*(simul_data[window, 0, 0:array_len]-centers[0,window])
	
	mean_pot[0, window] = np.mean(potentials[window, 0:array_len])

integ_param = 0
factor = 0
if fold_analyse != 'u':
	integ_param = constants[0]
	factor = 0.5
else:
	integ_param = centers[0]
	factor = -1

f=interp1d(integ_param, mean_pot[0],kind="cubic")

factor*intSimps(integ_param[0], integ_param[-1], 10000, 2)