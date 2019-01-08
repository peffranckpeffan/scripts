import numpy as np
from math import *
import os, errno, glob, shutil, json
import sys
import pprint as pp
from lib.py import util



# for i in range(9, 42):
# 	print i
# 	os.rename('../US/u'+str(i+3).zfill(2), '../US/u'+str(i).zfill(2))

a=0
for folder in os.listdir('../restraints/'):
	if os.path.isdir('../restraints/'+folder) and 'b' in folder:
		a = a + 1

print(a)




