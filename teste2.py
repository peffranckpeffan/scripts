import numpy as np
from math import *
import os, errno, glob, shutil, json
import sys
import pprint as pp
from lib.py import util

folders=["b", "l", "p", ""]
completed = 'false'
for i in range(0,15):
	file = open("../restraints/b"+str(i).zfill(2)+"/conf_rest-"+str(i).zfill(2)+".log")
	lines = file.readlines()
	for line in lines:
		if 'RITING VELOCITIES TO OUTPUT FILE AT STEP 1000000' in line:
			completed = 'true'
	print(completed + " "+ str(i).zfill(2))