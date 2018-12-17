import numpy as np
from math import *
import os, errno, glob, shutil, json

def copyAllFilesWith(source_dir, dest_dir, expression):
	files = glob.iglob(os.path.join(source_dir, expression))
	for file in files:
		print 'fasdfd'
		if (os.path.isfile(file)):
			shutil.copy2(file, dest_dir)

copyAllFilesWith('../common/','../equil/', '*eq*')

