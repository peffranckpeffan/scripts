import sys
import subprocess
import os
import time
from shutil import copyfile
import numpy as np
from math import *
import os, errno, glob, shutil, json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

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

z_positions = np.zeros([num_win, array_len], np.float64)
potentials = np.zeros([num_win, array_len], np.float64)
constants = np.zeros([1, num_win], np.float64)
z0 = np.zeros([1, num_win], np.float64)
average_pot = np.zeros([1, num_win], np.float64)

line_ref = 0

rest_info = []
if(fold_analyse == 'u'):
	rest_info = ['posit3', 2]
else:
	rest_info = ['rmsd2', 1]

for window in range(num_win):
	z_positions[window,0:array_len] = np.genfromtxt(dir_analyse+fold_analyse+str(window).zfill(2)+'/restraints.dat', usecols=(rest_info[1])).transpose()
	
	with open(dir_analyse+fold_analyse+str(window).zfill(2)+'/colvar.in') as file:
		lines = file.readlines()
		for line in lines:
			if rest_info[0] in line:
				line_ref = lines.index(line)
		z0[0, window] = float(lines[line_ref+2].split()[1])
		constants[0, window] = float(lines[line_ref+3].split()[1])

	if fold_analyse == 'u':	
		potentials[window, 0:array_len] = constants[0, window]*(z_positions[window,0:array_len]-z0[0,window])
	else:
		max_const = np.max(constants[0, window])
		potentials[window, 0:array_len] = max_const*(z_positions[window,0:array_len]-z0[0,window])**2
	
	average_pot[0, window] = np.mean(potentials[window, 0:array_len])

f=interp1d(constants[0]/max_const,average_pot[0],kind="cubic")

print (1/2)*intSimps(z0[0][0], z0[0][-1], 1000, 2)