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


# force_constant = 10
# z_windows = []
# k = []
# z0 = []
# line_num = 0
# u = []

# for i in range(0,47):
# 	z_windows.append(np.genfromtxt('../analysis/u'+str(i).zfill(2)+'/restraints.dat', usecols=(2)).transpose())

# 	with open('../analysis/u'+str(i).zfill(2)+'/colvar.in' % k, 'r') as file:
# 		lines = file.readlines()
# 		for line in lines:
# 			if 'posit3' in line:
# 				line_num = lines.index(line)
# 		k.append(float(lines[line_num+3].split()[1]))
#     	z0.append(float(lines[line_num+2].split()[1]))

# #print z_windows
# for i in range(0,47): 
# 	teste = []
# 	for z in z_windows[i]:
# 		teste.append(k[i]*(z-z0[i]))
# 	u.append(teste)

# # print u[0]
# # print len(u)
# def intSimps(a, b, n, step):
# 		h = (b - a) / n
# 		s = 0
# 		i = 0
# 		s = f(a) + f(b)

# 		for i in np.arange(1,n,step):
# 			s+= 4*f(a + i*h)
# 		for i in np.arange(2,n,step):
# 			s+= 2*f(a + i*h)
		
# 		s = (s*h)/3
# 		return s

# I=[]


# h = []
# for i in range(0,47):
# 	h.append(np.mean(u[i]))
# # plt.plot(z0, h)
# # plt.show()

# f=interp1d(z0,h,kind="cubic")

# file = open('forces.dat', 'w')
# # for i in range(len(z0)):
# # 	file.write("%15.7f %15.7f \n" % (z0[i], h[i]))
# print -1*intSimps(z0[0], z0[-1], 10000, 2)
# file.close()
# print np.trapz(z0,h)

# # for x in range(0,47):

# # 	f=interp1d(z0,h,kind="cubic")

# # 	I.append(-1*intSimps(z0[0], z0[-1], 1000, 2))

# # print np.sum(I)
# #print np.mean(I)

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
for window in range(num_win):
	z_positions[window,0:array_len] = np.genfromtxt(dir_analyse+fold_analyse+str(window).zfill(2)+'/restraints.dat', usecols=(2)).transpose()
	
	with open(dir_analyse+fold_analyse+str(window).zfill(2)+'/colvar.in') as file:
		lines = file.readlines()
		for line in lines:
			if 'posit3' in line:
				line_ref = lines.index(line)
		z0[0, window] = float(lines[line_ref+2].split()[1])
		constants[0, window] = float(lines[line_ref+3].split()[1])
    	
	potentials[window, 0:array_len] = (z_positions[window,0:array_len]-z0[0,window])

	average_pot[0, window] = constants[0, window]*np.mean(potentials[window, 0:array_len])

f=interp1d(z0[0],average_pot[0],kind="cubic")

print -1*intSimps(z0[0][0], z0[0][-1], 10000, 2)



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

simul_data = np.zeros([num_win, array_len], np.float64)
potentials = np.zeros([num_win, array_len], np.float64)
constants = np.zeros([1, num_win], np.float64)
centers = np.zeros([1, num_win], np.float64)
average_pot = np.zeros([1, num_win], np.float64)


rest_info = ['rmsd2', 2]
line_ref = 0
for window in range(num_win):
	simul_data[window,0:array_len] = np.genfromtxt(dir_analyse+fold_analyse+str(window).zfill(2)+'/restraints.dat', usecols=(rest_info[1])).transpose()

	with open(dir_analyse+fold_analyse+str(window).zfill(2)+'/colvar.in') as file:
		lines = file.readlines()
		for line in lines:
			if rest_info[0] in line:
				line_ref = lines.index(line)
		centers[0, window] = float(lines[line_ref+2].split()[1])
		constants[0, window] = float(lines[line_ref+3].split()[1])

	max_const = np.max(constants[0, window])
		
	potentials[window, 0:array_len] = (simul_data[window,0:array_len]-centers[0,window])**2
	
	average_pot[0, window] = np.mean(potentials[window, 0:array_len])

f=interp1d(constants[0], average_pot[0],kind="cubic")
# file = open('forces-r.dat', 'w')
# for i in range(num_win):
# 	file.write("%15.7f %15.7f \n" % (constants[0][i]/max_const, average_pot[0][i]))
# file.close()
print .5*intSimps(constants[0][0], constants[0][-1], 10000, 2)


# plt.plot(constants[0]/max_const, average_pot[0], 'o' )
# plt.show()