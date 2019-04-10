import sys
import subprocess
import os
import time
from shutil import copyfile
import numpy as np
from math import *
import os, errno, glob, shutil, json
from scipy.interpolate import interp1d


force_constant = 10
z_windows = []
k = []
z0 = []
line_num = 0
u = []

for i in range(0,47):
	z_windows.append(np.genfromtxt('../analysis/u'+str(i).zfill(2)+'/restraints.dat', usecols=(2)).transpose())

	with open('../analysis/u'+str(i).zfill(2)+'/colvar.in' % k, 'r') as file:
		lines = file.readlines()
		for line in lines:
			if 'posit3' in line:
				line_num = lines.index(line)
		k.append(float(lines[line_num+3].split()[1]))
    	z0.append(float(lines[line_num+2].split()[1]))

print z_windows
for i in range(0,47): 
	teste = []
	for z in z_windows[i]:
		teste.append(k[i]*(z-z0[i]))
	u.append(teste)

# print u[0]
# print len(u)
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

I=[]

for x in range(0,47):

	f=interp1d(z_windows[x],u[x],kind="cubic")


	I.append(-1*intSimps(z_windows[x][0], z_windows[x][len(z_windows[x])-1], 1000, 2))

print np.sum(I)
#print np.mean(I)
