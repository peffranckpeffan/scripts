import numpy as np
from math import *
import os, errno, glob, shutil, json
import sys
import pprint as pp
from collections import OrderedDict

# window_path  = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0]
# totrange = window_path[len(window_path)-1]-window_path[0]

# for i in range(7, 9):
# 	win_num = str(i).zfill(2)+'.1'
# 	files = os.listdir('../US/u'+win_num)
# 	for file in files:
# 		os.rename('../US/u'+win_num+'/'+file, '../US/u'+win_num+'/'+file.replace(win_num, str(i-1).zfill(2)+'.1')) 
# 	os.rename('../US/u'+win_num+'/', '../US/u'+str(i-1).zfill(2)+'.1/')

# dir_analyse = '../'+'US'+'/'
# teste=[]
# for window in os.listdir(dir_analyse):
# 	if os.path.isdir(dir_analyse+window) and 'u' in window:
# 		teste.append(window)
# #print teste
# #print(sorted(teste, key=lambda x: float(x[1:])))

# for index, window in enumerate(sorted(teste, key=lambda x: float(x[1:]))):
# 	print index, window.split('u')

# print len(window_path)
#util.call_subprocess("env totrange='"+str(totrange)+"' window_path='"+str(window_path).replace(' ', '').replace('[','').replace(']', '')+"' vmd -dispdev text -e ../scripts/lib-tcl/get-frames.tcl", "../US", True)

# teste = [[]]*47

# print teste

# s = [1, 2, 3]

# print s[-1]

# a = np.zeros([1,2])

# print a

# from pyexcel_ods import save_data

# data = OrderedDict()
# data.update({"Sheet 1": [
# 						['run', '0', '', 1, '', 2, '',3,'', 4], 
# 						['', 'MBAR', 'TI', 'MBAR', 'TI', 'MBAR', 'TI','MBAR', 'TI', 'MBAR', 'TI'], 
# 						['Protein Conformational Bulk', str(energies['b'][0]) + '+-' str(energies['b'][0]), 6],
# 						['Protein Conformational Site', str(energies['p'][0]) + '+-' str(energies['p'][0]), 6],
# 						['Ligand Conformational Bulk', str(energies['l'][0]) + '+-' str(energies['l'][0]), 6],
# 						['Ligand Conformational Site', str(energies['r'][0]) + '+-' str(energies['r'][0]), 6],
# 						['Ligand Orientational Bulk', str(energies['b'][0]) + '+-' str(energies['b'][0]), 6],
# 						['Ligand Orientational Site', str(energies['o'][0]) + '+-' str(energies['o'][0]), 6],
# 						['Ligand Translational Bulk', str(energies['b'][0]) + '+-' str(energies['b'][0]), 6],
# 						['Ligand Translational Site', str(energies['t'][0]) + '+-' str(energies['t'][0]), 6],
# 						['Potential of Mean Force', str(energies['u'][0]) + '+-' str(energies['u'][0]), 6]
# 						['Total', 5, 6]
# 	]})

# save_data("your_file.ods", data)

run = 0

x = y = z = np.arange(0.0,5.0,1.0)
teste = np.zeros([10,10])
teste = np.genfromtxt('test.out', usecols=(0,1,2,3,4,5,6,7,8,9)).transpose()


for i in range():
	teste[run+1,0:10] = 
	teste[run+2,0:10] = 

np.savetxt('test.out', )   # x,y,z equal sized 1D arrays