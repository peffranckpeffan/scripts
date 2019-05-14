#!/usr/bin/python

import sys
import subprocess
import os
import time
from shutil import copyfile
from pathlib import Path
import numpy as np
from math import *
import os, errno, glob, shutil, json
from pprint import pprint
from lib.py import util
from lib.py import analysis
from lib.py import integration_methods as im

#equil - python index.py equil {exec}
#smd   - python index.py smd {exec}
#restraints - python index.py restraints [folder to be executed (t, o, r, p, b, all)] [phase] {exec}
#us - python index.py us [stage] {exec}
#analysis - python index.py analysis [(us or restraints)]->[specific restraint?] [phase]

#NEEDS TO BE INSTALLED
#akima, pymbar, scipy, timeseries, configparser, pathlib
goFoward=True
stage=""
res_folder=""
ifExec = 0

stages=['equil', 'smd', 'restraints', 'us', 'analysis', 'init']
folders_rest=['b', 'l', 'o', 'p', 'r', 't']

try:
	stage=sys.argv[1]
except:
	sys.exit("It's necessary to inform de stage to be generated.")

if stage not in stages:
	sys.exit('Stage of the analysis informed don\'t exist.')

else:

	if (stage == 'equil'):
		if ("exec" in sys.argv):
			ifExec = 1

	if (stage == 'smd'):
		if ("exec" in sys.argv):
			ifExec = 1
	
	if (stage == 'restraints'):
		if ("exec" in sys.argv and len(sys.argv) < 5) or ("exec" not in sys.argv and len(sys.argv) < 4):
			sys.exit("Missing options. Check sintax")

		try:
			res_folder=sys.argv[2]
		except:
			sys.exit("It's necessary to inform the restraint to be generated.")

		try:
			phase=sys.argv[3]
		except:
			sys.exit("It's necessary to inform de phase of the simulation to be generated.")

		if (res_folder not in folders_rest and res_folder != 'all' ):
			sys.exit("Folder especificied do not exist. Please inform: t, o, r, p, b, all")
		if ("exec" in sys.argv):                                                                                                        
		 	ifExec = 1

	if (stage == 'us'):
		if ("exec" in sys.argv and len(sys.argv) < 4) or ("exec" not in sys.argv and len(sys.argv) < 2):
			sys.exit("Missing options. Check sintax.")

		try:
			phase=sys.argv[2]
		except:
			sys.exit("It's necessary to inform de phase of the simulation to be generated.")

		if ("exec" in sys.argv):
			ifExec = 1
	
	if (stage == 'analysis'):
		try:
			stageAnalysis=sys.argv[2]
		except:
			sys.exit("It's necessary to inform de stage to be analised.")

		if (stageAnalysis != 'restraints' and stageAnalysis != 'us'):
			sys.exit('Stage to be analized is not valid.')

		if (stageAnalysis=='restraints' and len(sys.argv) < 5) or (stageAnalysis=='us' and len(sys.argv) < 4):
			sys.exit("Missing options. Check sintax.")
		
		phase_position=3
		if (stageAnalysis == 'restraints'):
			try:
				res_folder=sys.argv[3]
			except:
				sys.exit("It's necessary to inform the folder to be analized. Please inform: t, o, r, p, b, all")
			
			phase_position=4

			if (res_folder not in folders_rest and res_folder != 'all' ):
				sys.exit("Folder especificied do not exist. Please inform: t, o, r, p, b, all")			

		try:
			phase=sys.argv[phase_position]
		except:
			sys.exit("It's necessary to inform de phase (0, 1, 2, ...) of the simulation to be analised.")


if (stage == "init"):
	createDir("../common")
	copyAllFilesWith('lib/param','../common', '*')
	createDir("../equil")
	copyAllFilesWith('../common','../equil', '*eq*')
	createDir("../restraints")
	copyAllFilesWith('../common','../restraints', '*rest*')
	createDir("../US")
	copyAllFilesWith('../common','../US', '*us*')

elif (stage == "equil"):
	util.createDir("../equil")
	util.copyAllFilesWith('lib/colv-files','../equil', '*eq*')
	util.copyAllFilesWith('lib/conf-files','../equil', '*eq*')

	#Mount the system
	util.call_subprocess("vmd -dispdev text -e ../scripts/lib/tcl/setup-equil.tcl", "../equil", True)

	#Update the dummy atom
	util.update_dummy('equil', ['colv-equil'], "../common/system.pdb", "segname A and backbone", stage)

	copyAllFilesWith('../common','../equil', '*')

	if (ifExec==1):
		file_eq=Path('../equil/out_eq.dcd')
		if (file_eq.exists()):
			util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq2 > conf-eq2.log", "../equil", True)
		else:
			util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq > conf-eq.log", "../equil", True)
			util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq2 > conf-eq2.log", "../equil", True)
		

elif (stage == "smd"):
	util.createDir("../SMD")
	util.copyAllFilesWith('lib/colv-files','../SMD', '*smd*')
	util.copyAllFilesWith('lib/conf-files','../SMD', '*smd*')

	#Mount the system
	util.call_subprocess("vmd -dispdev text -e ../scripts/lib/tcl/setup-smd.tcl", "../SMD", True)

	#Update dummyatom
	util.update_dummy('SMD', ['colv-smd'], "../equil/out_eq2.restart.coor", "segname B and backbone", stage)

	util.update_dummy('SMD', ['colv-smd'], "../equil/out_eq2.restart.coor", "segname A and backbone", stage)

	copyAllFilesWith('../common','../SMD', '*')
	copyAllFilesWith('../equil','../SMD', '*eq2_restart.coor*')

	if (ifExec==1 and Path('../equil/out_eq2.dcd').exists()):
		util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-smd > conf-smd.log", "../SMD", True)

elif (stage == "us"):
	if not(Path('../US').exists()):
		util.createDir("../US")
		util.copyAllFilesWith('lib/colv-files','../US', '*us*')
		util.copyAllFilesWith('lib/conf-files','../US', '*us*')

	window_path  = [4.0, 4.5, 5.0, 5.25, 5.5, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0]
	totrange = window_path[len(window_path)-1]-window_path[0]
	
	util.call_subprocess("vmd -dispdev text -e ../scripts/lib/tcl/setup-smd.tcl", "../US", True)

	util.update_dummy('US', ['basic_colv-us'], "../common/system.pdb", "segname A and backbone", 'equil')
	util.update_dummy('US', ['basic_colv-us'], "../equil/out_eq2.restart.coor", "segname B and backbone", 'smd')

	if(not(Path('../US/u40/out_min-40.restart.coor').exists())):
		#PREPARING MINIMIZATION
		steps=500000
		#initial window in posit_wind + space_bet_wind (3.5)
		i = 0
		for position in window_path:

			window=str(i).zfill(2)
			
			inputt="min-"+str(window)
			input_pr="smd-"+str(window)
			conf_name="conf_"+inputt
			colv_inp="colv-"+str(window)

			print("conf_name: " + conf_name + "; input: " + "; " + inputt +"," + input_pr + "; window: " + str(window))

			util.createDir("../US/u"+str(window))

			conf_file = open("../US/u"+str(window)+"/"+conf_name, 'w')
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "# INPUT AND OUTPUT FILES                           ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "set input " + inputt + "\n" )
			conf_file.write( "set input_pr " + input_pr + "\n" )
			conf_file.write( "set colv_inp " + colv_inp + "\n" )
			conf_file.write( "set inputname   ./out_$input_pr" + "\n" )
			conf_file.write( "set outputname  ./out_$input" + "\n" )
			conf_file.write( "bincoordinates  $inputname.restart.coor" + "\n" )
			conf_file.write( "#binvelocities   $inputname.restart.vel" + "\n" )
			conf_file.write( "extendedSystem  ../../SMD/out_smd.restart.xsc" + "\n" )
			conf_file.write( "set ref_umb     ../refumb0.pdb" + "\n" )
			conf_file.write( "coordinates     ../../common/system.pdb" + "\n" )
			conf_file.write( "structure       ../../common/system.psf" + "\n" )
			conf_file.write( "" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "## INPUT SETTINGS                                   ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "set temperature    300" + "\n" )
			conf_file.write( "set cons  1" + "\n" )
			conf_file.write( "set min   1" + "\n" )
			conf_file.write( "set pres   1" + "\n" )
			conf_file.write("temperature  $temperature" + "\n")
			conf_file.write( "" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "## JOB DESCRIPTION                                  ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "" )
			util.copyFileLines("lib/conf-files/basic_conf-us", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();

			tcl_file = open("../US/u"+str(window)+"/"+colv_inp, "w")
			util.copyFileLines("lib/colv-files/basic_colv-us", tcl_file)
			tcl_file.write( "" )
			tcl_file.write( "harmonic { " + "\n")
			tcl_file.write( "  name posit3 \n" )
			tcl_file.write( "  colvars posit-Z-lig \n" )
			tcl_file.write( "  centers " +  str(position) + "\n" )
			tcl_file.write( "  forceConstant 10.0 \n" )
			tcl_file.write( "} " + "\n")
			tcl_file.write( "" )
			tcl_file.close();
			i=i+1

	if(not(Path('../US/u40/out_smd-40.restart.coor').exists())):
		util.call_subprocess("env totrange='"+str(totrange)+"' window_path='"+str(window_path).replace(' ', '').replace('[','').replace(']', '')+"' vmd -dispdev text -e ../scripts/lib/tcl/get-frames.tcl", "../US", True)



	#PREPARING TO RUN US
	steps=1000000
	#initial window in posit_wind + space_bet_wind (3.5)
	i = 0

	for position in window_path:
		window=str(i).zfill(2)

		if phase == '0':
			inputt="run-"+str(window)
			input_pr="min-"+str(window)
		elif phase == '1':
			inputt="run1-"+str(window)
			input_pr="run-"+str(window)
		else:
			inputt="run"+phase+"-"+str(window)
			input_pr="run"+str(int(phase)-1)+"-"+str(window)
		
		conf_name="conf_"+inputt
		colv_inp="colv-"+str(window)

		print("conf_name: " + conf_name + "; input: " + "; " + inputt +"," + input_pr + "; window: " + str(window))

		if(not(Path('../US/u'+str(window)+'/out_smd-'+str(window)+'.restart.coor').exists())):
			shutil.move('../US/out_smd-'+str(window)+'.restart.coor', "../US/u"+str(window))

		copyAllFilesWith('../common','../restraints/'+str(window), '*')
		copyAllFilesWith('../SMD/out_smd.restart.xsc"','../restraints/'+str(window), '*')
		copyAllFilesWith('../restraints/refumb0.pdb"','../restraints/'+str(window), '*')
		copyAllFilesWith('../restraints/atoms.pdb"','../restraints/'+str(window), '*')

		conf_file = open("../US/u"+str(window)+"/"+conf_name, 'w')
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "# INPUT AND OUTPUT FILES                           ##" + "\n" )
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "set input " + inputt + "\n" )
		conf_file.write( "set input_pr " + input_pr + "\n" )
		conf_file.write( "set colv_inp " + colv_inp + "\n" )
		conf_file.write( "set inputname   ./out_$input_pr" + "\n" )
		conf_file.write( "set outputname  ./out_$input" + "\n" )
		conf_file.write( "bincoordinates  $inputname.restart.coor" + "\n" )
		conf_file.write( "#binvelocities   $inputname.restart.vel" + "\n" )
		conf_file.write( "extendedSystem  ./SMD/out_smd.restart.xsc" + "\n" )
		conf_file.write( "set ref_umb     ./refumb0.pdb" + "\n" )
		conf_file.write( "coordinates     ./system.pdb" + "\n" )
		conf_file.write( "structure       ./system.psf" + "\n" )
		conf_file.write( "" )
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "## INPUT SETTINGS                                   ##" + "\n" )
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "set temperature    300" + "\n" )
		conf_file.write( "set cons  1" + "\n" )
		conf_file.write( "set min   0" + "\n" )
		conf_file.write( "set pres   1" + "\n" )
		conf_file.write("temperature  $temperature" + "\n")
		conf_file.write( "" )
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "## JOB DESCRIPTION                                  ##" + "\n" )
		conf_file.write( "######################################################" + "\n" )
		conf_file.write( "" )
		util.copyFileLines("lib/conf-files/basic_conf-us", conf_file)
		conf_file.write("run "+ str(steps))
		conf_file.close();
		i=i+1

	print("PHASE: "+ phase)
	if ifExec==1:
		if (phase=='0'):
			process=["min", "run"]
			print("The minization of the system will be perfomed followed by the run "+phase+" of the PMF.")
		elif(phase!='0' and not(Path('../US/u40/out_min-40.restart.coor').exists())):
			process=["min", "run"]
			print("The minization of the system will be perfomed followed by the run 0 of the PMF.")
		elif(phase =='1' and Path('../US/u40/out_run-40.restart.coor').exists()):
			process=["run1"]
		elif(phase!='0' and Path('../US/u40/out_run'+str(int(phase)-1)+'-40.restart.coor').exists()):
			process=["run"+str(phase)]
		else:
			print("The pointed phase of the PMF simulations does not fit with the previous runs, check yours files and try again.")
			quit()
		for x in range(len(process)):
			for i in range(0,len(window_path)):
				window=str(i).zfill(2)
				print(process[x] + " " + window)
				util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_"+process[x]+"-"+str(window)+" > conf_"+process[x]+"-"+str(window)+".log", "../US/u"+str(window), True)

#CALCULATE RESTRICTIONS
elif (stage == 'restraints'):
	if not(Path('../restraints').exists()):
		util.createDir("../restraints")
		util.copyAllFilesWith('lib/colv-files','../restraints', '*basic_colv*')
		util.copyAllFilesWith('lib/conf-files','../restraints', '*basic_conf-rest*')

	steps=1000000
	restraints=json.load(open("restraints.json"))["restraints"]
	
	util.call_subprocess("vmd -dispdev text -e ../scripts/lib/tcl/setup-smd.tcl", "../restraints", True)

	basic_colv_list=['basic_colv-rest-bulk', 'basic_colv-rest-lrmsd', 'basic_colv-rest-orient', 'basic_colv-rest-prmsd', 'basic_colv-rest-trans']
	util.update_dummy('restraints', basic_colv_list, "../common/system.pdb", "segname A and backbone", 'equil')
	util.update_dummy('restraints', basic_colv_list, "../equil/out_eq2.restart.coor", "segname B and backbone", 'smd')

	multipliers=[0.00, 0.40, 0.80, 1.60, 2.40, 4.00, 5.50, 8.65, 11.80, 18.10, 24.40, 37.00, 49.60, 74.80, 100.00]
	j="0"

	count_rest=0
	while  len(restraints) > count_rest:
		#ft=float(restraints[count_rest]["forceConstant"])

		counter = 0
		while  len(multipliers) > counter:

			ct=multipliers[counter]
			#fw=ft*ct/100

			if (counter < 10):
				nr=j+str(counter)
			else:
				nr=counter


			if phase=='0':
				if restraints[count_rest]["folder"] == 'b':
					inputt='smd'
					input_dir='../../SMD'
				else:
					inputt="eq2"
					input_dir='../../equil'
				output="rest-" + str(nr)
			elif phase=='1':
				inputt="rest-" + str(nr)
				input_dir='.'
				output="rest1-" + str(nr)
			else:
				inputt="rest"+str(int(phase)-1)+"-" + str(nr)
				input_dir='.'
				output="rest"+phase+"-" + str(nr)	


			# if phase == "0":
			# 	phase = ""

			conf_name="conf_" + output
			colv_inp="colv-" + str(nr)
			
			if not(Path("../restraints/"+restraints[count_rest]["folder"] + str(nr)).exists()):
				util.createDir("../restraints/"+restraints[count_rest]["folder"] + str(nr))

			print("conf_name: " + conf_name + "; force_const: " + str(restraints[count_rest]["properties"][0]['forceConstant']*ct/100) +"; " + restraints[count_rest]["folder"]  + str(nr))

			conf_file = open("../restraints/"+restraints[count_rest]["folder"]+str(nr)+"/"+conf_name, 'w')
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "# INPUT AND OUTPUT FILES                           ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "set input " + inputt + "\n" )
			conf_file.write( "set output " + output + "\n" )
			conf_file.write( "set colv_inp " + colv_inp + "\n" )
			conf_file.write( "set inputname   "+input_dir+"/out_$input" + "\n" )
			conf_file.write( "set outputname  ./out_$output" + "\n" )
			conf_file.write( "bincoordinates  $inputname.restart.coor" + "\n" )
			conf_file.write( "binvelocities   $inputname.restart.vel" + "\n" )
			conf_file.write( "extendedSystem  $inputname.restart.xsc" + "\n" )
			conf_file.write( "set ref_umb     ./refumb0.pdb" + "\n" )
			conf_file.write( "coordinates     ./system.pdb" + "\n" )
			conf_file.write( "structure       ./system.psf" + "\n" )
			conf_file.write( "" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "## INPUT SETTINGS                                   ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "set temperature    300" + "\n" )
			conf_file.write( "set cons  1" + "\n" )
			conf_file.write( "set min   0" + "\n" )
			conf_file.write( "set pres   1" + "\n" )
			conf_file.write( "" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "## JOB DESCRIPTION                                  ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "" )
			util.copyFileLines("lib/conf-files/basic_conf-rest", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();


			colv_file = open("../restraints/"+restraints[count_rest]["folder"]+str(nr)+"/"+colv_inp, "w")
			util.copyFileLines("lib/colv-files/basic_colv-rest-" + restraints[count_rest]["type"], colv_file)
			counter_rest2=0
			
			while len(restraints[count_rest]["properties"]) > counter_rest2:
				colv_file.write( "" )
				colv_file.write( "harmonic { " + "\n")
				colv_file.write( "  name " + restraints[count_rest]["properties"][counter_rest2]["name"] + "\n" )
				colv_file.write( "  colvars " + restraints[count_rest]["properties"][counter_rest2]["colvars"] + "\n" )
				colv_file.write( "  centers " +  restraints[count_rest]["properties"][counter_rest2]["centers"] + "\n" )
				colv_file.write( "  forceConstant " + str(float(restraints[count_rest]["properties"][counter_rest2]["forceConstant"])*ct/100) + "\n" )
				colv_file.write( "} " + "\n")
				colv_file.write( "" )

				counter_rest2=counter_rest2+1

			colv_file.close();


			counter=counter+1
		count_rest=count_rest+1

	if ifExec==1:

		if phase == "0":
			phase = ""

		folders=["b", "o", "p", "t", "r"]
		#folders=["t", "o", "r"]
		if res_folder not in folders and res_folder != "all":
			print("Specified folder does not exist.")
		else:
			if res_folder in folders:
				folders=[res_folder]

			for x in range(len(folders)):
				for i in range(0,15):
					window=str(i).zfill(2)
					util.call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_rest"+phase+"-"+str(window)+" > conf_rest"+phase+"-"+str(window)+".log", "../restraints/"+folders[x]+str(window), True)


elif stage=="analysis":
	if (phase == '0'):
		phase = ""

	if (stageAnalysis=="us"):
		prefix = "run" + phase
		folders = ["u"]
		stageAnalysis = stageAnalysis.upper()
	else:
		prefix ="rest" + phase
		if (res_folder == 'all'):
			folders = folders_rest
		else:
			folders = [res_folder]

	dir_analyse = '../'+stageAnalysis+'/'
	windows_analyse = []
	for prefix_folder in folders:
		folder_to_count = prefix_folder

		if (prefix_folder == 'l'):
			folder_to_count = 'b'

		for window in os.listdir(dir_analyse):
			
			if os.path.isdir(dir_analyse+window) and folder_to_count in window:
				windows_analyse.append(window)

		windows_analyse = sorted(windows_analyse, key=lambda x: float(x[1:]))
		for index, window in enumerate(windows_analyse):
			
			folder_from = window
	 		folder_for = prefix_folder + str(index).zfill(2)

			util.createDir("../analysis/"+folder_for)
			util.copyfile("../"+stageAnalysis+"/"+folder_from+"/out_"+prefix+"-"+window.split(folder_to_count)[1]+".colvars.traj", "../analysis/"+folder_for+"/restraints.dat")
			util.copyfile("../"+stageAnalysis+"/"+folder_from+"/colv-"+window.split(folder_to_count)[1], "../analysis/"+folder_for+"/colvar.in")

		print("Analyzing folder " + prefix_folder)
		util.call_subprocess("python ../scripts/FE-MBAR.py "+prefix_folder+" 298 "+ phase, "../analysis", True)
		util.call_subprocess("python ../scripts/FE-MBAR-ns.py "+prefix_folder+" 298 "+ phase, "../analysis", True)

		del windows_analyse[:]

	energies = []
	thermo_intg = []
	fileType = ['allsp', 'subs']
	total = []
	if(Path("../analysis/"+fileType[0]+"-t.100.dat").exists() and Path("../analysis/"+fileType[0]+"-u.100.dat").exists()):
		for folder in ['b', 'l', 'o', 'p', 'r', 't', 'u']:
			vals = []
			for y in range(1,3):#valor absoluto posicao 1 e seu erro 2
				vals.append(util.getEnergy(fileType[y-1]+"-"+folder+".100.dat", y))
			vals.append(im.thermodynamic_integration(folder))
			energies.append(vals)

		energies.insert(4, [6.72, .0, 6.72])
		energies.insert(6, [5.27, .0, 5.27])

		total_mbar = energies[0][0] - energies[3][0] + energies[1][0] - energies[5][0] + energies[4][0]  - energies[2][0] + energies[6][0]  - energies[7][0] - energies[8][0]
		total_ti = energies[0][2] - energies[3][2] + energies[1][2] - energies[5][2] + energies[4][2]  - energies[2][2] + energies[6][2]  - energies[7][2] - energies[8][2]
		
		energies.append([total_mbar, .0, total_ti])

		

		# teste = ['Protein Conformational Bulk', 'Protein Conformational Site', 'Ligand Conformational Bulk', 'Ligand Conformational Site', 'Ligand Orientational Bulk', 'Ligand Orientational Site', 'Ligand Translational Bulk', 'Ligand Translational Site', 'PMF', 'Total']

		# columns[0:10,0] = teste[0:10]
		if not(Path('../analysis/RESULT.out').exists()):
			columns = np.zeros([10,15])
			np.savetxt('../analysis/RESULT.out', (columns[0:10,0:15]), delimiter='\t', fmt='%10.3f ')

		data = np.genfromtxt('../analysis/RESULT.out', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
		
		if phase == '':
			phase = '0'

		x = int(phase)
		for i in range(0,10):
			data[i][x*3:(x+1)*3] = energies[i][0:3]

		np.savetxt('../analysis/RESULT.out', (data[0:10,0:15]), delimiter='\t', fmt='%10.3f ')