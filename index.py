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

#equil - python index.py equil {exec}
#smd   - python index.py smd {exec}
#restraints - python index.py restraints [folder to be executed (t, o, r, p, b, all)] [phase] {exec}
#us - python index.py us [stage] {exec}
#analysis - python index.py analysis [(us or restraints)]->[specific restraint?] [phase]
goFoward=True
stage=""
res_folder=""

stages=['equil', 'smd', 'restraints', 'us', 'analysis', 'init']
folders_rest=['b', 'l', 'o', 'p', 'r', 't']

try:
	stage=sys.argv[1]
except:
	sys.exit("It's necessary to inform de stage to be generated.")

if stage not in stages:
	sys.exit('Stage of the analysis informed don\'t exist.')

else:
	
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

	if (stage == 'us'):
		if ("exec" in sys.argv and len(sys.argv) < 4) or ("exec" not in sys.argv and len(sys.argv) < 3):
			sys.exit("Missing options. Check sintax.")

		try:
			phase=sys.argv[2]
		except:
			sys.exit("It's necessary to inform de phase of the simulation to be generated.")
	
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
	createDir("../equil")
	copyAllFilesWith('../common','../equil', '*eq*')
	createDir("../restraints")
	copyAllFilesWith('../common','../restraints', '*rest*')
	createDir("../SMD")
	copyAllFilesWith('../common','../SMD', '*smd*')
	createDir("../US")
	copyAllFilesWith('../common','../US', '*us*')

elif (stage == "equil"):

	#Mount the system
	call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-equil.tcl", "../equil", True)

	#Update the dummy atom
	update_dummy('equil', ['colv-equil'], "../common/solvate.pdb", "segname A and backbone", stage)

	if (ifExec=="exec"):
		file_eq=Path('../equil/out_eq.dcd')
		if (file_eq.exists()):
			call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq2 > conf-eq2.log", "../equil", True)
		else:
			call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq > conf-eq.log", "../equil", True)
			call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-eq2 > conf-eq2.log", "../equil", True)
		

elif (stage == "smd"):

	#Mount the system
	call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../SMD", True)

	#Update dummyatom
	update_dummy('SMD', ['colv-smd.inp'], "../equil/out_eq2.restart.coor", "segname B and backbone", stage)

	if (ifExec=="exec" and Path('../equil/out_eq.dcd').exists()):
		call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf-smd > conf-smd.log", "../SMD", True)

elif (stage == "us"):

	call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../US", True)

	update_dummy('common', ['basic_colv-us'], "../common/solvate.pdb", "segname A and backbone", 'equil')
	update_dummy('common', ['basic_colv-us'], "../equil/out_eq2.restart.coor", "segname B and backbone", 'smd')

	if(not(Path('../US/u41/out_min-41.restart.coor').exists())):
		#PREPARING MINIMIZATION
		steps=500000
		#initial window in posit_wind + space_bet_wind (3.5)
		posit_wind=3.0
		space_bet_wind=0.5

		for i in range(42):
			window=str(i).zfill(2)
			
			inputt="min-"+str(window)
			input_pr="smd-"+str(window)
			conf_name="conf_"+inputt
			colv_inp="colv-"+str(window)
			posit_wind = posit_wind + space_bet_wind

			print("conf_name: " + conf_name + "; input: " + "; " + inputt +"," + input_pr + "; window: " + str(window))

			createDir("../US/u"+str(window))

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
			conf_file.write( "coordinates     ../../common/solvate.pdb" + "\n" )
			conf_file.write( "structure       ../../common/solvate.psf" + "\n" )
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
			copyFileLines("../common/basic_conf-us", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();

			tcl_file = open("../US/u"+str(window)+"/"+colv_inp, "w")
			copyFileLines("../common/basic_colv-us", tcl_file)
			tcl_file.write( "" )
			tcl_file.write( "harmonic { " + "\n")
			tcl_file.write( "  name posit3 \n" )
			tcl_file.write( "  colvars posit-Z-lig \n" )
			tcl_file.write( "  centers " +  str(posit_wind) + "\n" )
			tcl_file.write( "  forceConstant 10.0 \n" )
			tcl_file.write( "} " + "\n")
			tcl_file.write( "" )
			tcl_file.close();

	if(not(Path('../US/u41/out_smd-41.restart.coor').exists())):
		call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/get-frames.tcl", "../US", True)


	#PREPARING TO RUN US
	steps=1000000
	#initial window in posit_wind + space_bet_wind (3.5)
	posit_wind=3.0
	space_bet_wind=0.5

	for i in range(0,42):
		window=str(i).zfill(2)

		inputt="run-"+str(window)
		input_pr="min-"+str(window)
		conf_name="conf_"+inputt
		colv_inp="colv-"+str(window)
		posit_wind = posit_wind + space_bet_wind

		print("conf_name: " + conf_name + "; input: " + "; " + inputt +"," + input_pr + "; window: " + str(window))

		if(not(Path('../US/u'+str(window)+'/out_smd-'+str(window)+'.restart.coor').exists())):
			shutil.move('../US/out_smd-'+str(window)+'.restart.coor', "../US/u"+str(window))

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
		conf_file.write( "coordinates     ../../common/solvate.pdb" + "\n" )
		conf_file.write( "structure       ../../common/solvate.psf" + "\n" )
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
		copyFileLines("../common/basic_conf-us", conf_file)
		conf_file.write("run "+ str(steps))
		conf_file.close();

	print("PHASE: "+ phaseUs)
	if ifExec=="exec":
		if (phaseUs=='0'):
			process=["min", "run"]
			print("The minization of the system will be perfomed followed by the run "+phaseUs+" of the PMF.")
		elif(phaseUs!='0' and not(Path('../US/u41/out_min-41.restart.coor').exists())):
			process=["min", "run"]
			print("The minization of the system will be perfomed followed by the run 0 of the PMF.")
		elif(phaseUs!='0' and Path('../US/u41/out_run'+str(int(phaseUs)-1)+'-41.restart.coor').exists()):
			process=["run"+str(phaseUs)]
		else:
			print("The pointed phase of the PMF simulations does not fit with the previous runs, check yours files and try again.")
			quit()
		for x in range(len(process)):
			for i in range(0,42):
				window=str(i).zfill(2)
				print(process[x])
				call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_"+process[x]+"-"+str(window)+" > conf_"+process[x]+"-"+str(window)+".log", "../US/u"+str(window), True)

#CALCULATE RESTRICTIONS
elif (stage == 'restraints'):

	steps=1000000
	restraints=json.load(open("restraints.json"))["restraints"]
	
	call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../restraints", True)

	basic_colv_list=['basic_colv-rest-bulk', 'basic_colv-rest-lrmsd', 'basic_colv-rest-orient', 'basic_colv-rest-prmsd', 'basic_colv-rest-trans']
	update_dummy('restraints', basic_colv_list, "../common/solvate.pdb", "segname A and backbone", 'equil')
	update_dummy('restraints', basic_colv_list, "../equil/out_eq2.restart.coor", "segname B and backbone", 'smd')

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

			pr="eq2"
			inputt="rest-" + str(nr)
			input_pr=pr
			#input_fr="rest-" + str(fr)
			conf_name="conf_" + str(inputt)
			colv_inp="colv-" + str(nr)

			createDir("../restraints/"+restraints[count_rest]["folder"] + str(nr))

			print("conf_name: " + conf_name + "; force_const: " + str(restraints[count_rest]["properties"][0]['forceConstant']) +"; " + restraints[count_rest]["folder"]  + str(nr))

			conf_file = open("../restraints/"+restraints[count_rest]["folder"]+str(nr)+"/"+conf_name, 'w')
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "# INPUT AND OUTPUT FILES                           ##" + "\n" )
			conf_file.write( "######################################################" + "\n" )
			conf_file.write( "set input " + inputt + "\n" )
			conf_file.write( "set input_pr " + input_pr + "\n" )
			conf_file.write( "set colv_inp " + colv_inp + "\n" )
			conf_file.write( "set inputname   ../../equil/out_$input_pr" + "\n" )
			conf_file.write( "set outputname  ./out_$input" + "\n" )
			conf_file.write( "bincoordinates  $inputname.restart.coor" + "\n" )
			conf_file.write( "binvelocities   $inputname.restart.vel" + "\n" )
			conf_file.write( "extendedSystem  $inputname.restart.xsc" + "\n" )
			conf_file.write( "set ref_umb     ../../SMD/refumb0.pdb" + "\n" )
			conf_file.write( "coordinates     ../../common/solvate.pdb" + "\n" )
			conf_file.write( "structure       ../../common/solvate.psf" + "\n" )
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
			copyFileLines("../common/basic_conf-rest", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();


			colv_file = open("../restraints/"+restraints[count_rest]["folder"]+str(nr)+"/"+colv_inp, "w")
			copyFileLines("../common/basic_colv-rest-" + restraints[count_rest]["type"], colv_file)
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

	if ifExec=="exec":

		folders=["t", "o", "r", "p", "b"]
		if res_folder not in folders and res_folder != "all":
			print("Specified folder does not exist.")
		else:
			if res_folder in folders:
				folders=[res_folder]

			for x in range(len(folders)):
				for i in range(0,15):
					window=str(i).zfill(2)
					call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_rest-"+str(window)+" > conf_rest-"+str(window)+".log", "../restraints/"+folders[x]+str(window), True)


elif stage=="analysis":
	if (phase == '0'):
		phase = ""

	if (stageAnalysis=="us"):
		prefix = "run" + phase
		folders = ["u"]
		stageAnalysis = stageAnalysis.upper()
		init_fold=1
	else:
		prefix ="rest" + phase
		init_fold=0
		if (res_folder == 'all'):
			folders = folders_rest
		else:
			folders = [res_folder]


	for x in range(len(folders)):
		num_folder = 0
		folder_to_count = folders[x]
		
		if (folders[x] == 'l'):
			folder_to_count = 'b'
		
		while os.path.exists('../'+stageAnalysis+'/'+folder_to_count+str(num_folder).zfill(2)):
			num_folder = num_folder + 1
		
		for i in range(init_fold, num_folder):
			window=str(i).zfill(2)

			folder_copy = folders[x]+window
			folder_for = folder_copy

			if (folders[x] == 'l'):
				folder_copy = 'b' + window

			elif (folders[x] == 'u'):
				folder_for = folders[x]+str(i-1).zfill(2)

			util.createDir("../analysis/"+folder_for)

			util.copyfile("../"+stageAnalysis+"/"+folder_copy+"/out_"+prefix+"-"+window+".colvars.traj", "../analysis/"+folder_for+"/restraints.dat")
			util.copyfile("../"+stageAnalysis+"/"+folder_copy+"/colv-"+window, "../analysis/"+folder_for+"/colvar.in")

		print("Analyzing folder " + folders[x])
		util.call_subprocess("python ../scripts/FE-MBAR.py "+folders[x]+" 298 "+phase, "../analysis", True)
		util.call_subprocess("python ../scripts/FE-MBAR-ns.py "+folders[x]+" 298 "+phase, "../analysis", True)

	#writing relevant data to a file
	K_bT=((1.381e-23 * 6.022e23) / (4.184 * 1000.0))*298
	fileType = ['allsp', 'subs']
	folders = folders_rest.append('u')
	
	if(Path("../analysis/"+fileType[0]+"-t.100.dat").exists() and Path("../analysis/"+fileType[0]+"-u.100.dat").exists()):	
		freeEnergies = open("../analysis/RESULT.dat", "w")

		total = 0
		energies = {}
		
		for folder in ['b', 'l', 'o', 'p', 'r', 't', 'u']:
			vals = []
			for y in range(1,3):#valor absoluto posicao 1 e seu erro 2
				vals.append(util.getEnergy(fileType[0]+"-"+folder+".100.dat", y))
			energies[folder] = vals
			 
		freeEnergies.write("Protein Conformational Bulk	%10.5f +- %10.5f \n" % (energies['b'][0], energies['b'][1])) 
		freeEnergies.write("Protein Conformational Site	%10.5f +- %10.5f \n" % (energies['p'][0], energies['p'][1]))	
		freeEnergies.write("Ligand Conformational Bulk	%10.5f +- %10.5f \n" % (energies['l'][0], energies['l'][1])) 
		freeEnergies.write("Ligand Conformational Site	%10.5f +- %10.5f \n" % (energies['r'][0], energies['r'][1])) 
		freeEnergies.write("Ligand Orientational Bulk		6.72 \n")
		freeEnergies.write("Ligand Orientational Site	%10.5f +- %10.5f \n" % (energies['o'][0], energies['o'][1])) 
		freeEnergies.write("Ligand Translational Bulk		5.27 \n")
		freeEnergies.write("Ligand Translational Site	%10.5f +- %10.5f \n" % (energies['t'][0], energies['t'][1])) 
		freeEnergies.write("PMF							%10.5f +- %10.5f \n\n" % (energies['u'][0], energies['u'][1])) 
		
		total = energies['b'][0] - energies['p'][0] + energies['l'][0] - energies['r'][0] + 6.72 - energies['o'][0] + 5.27 - energies['t'][0] - energies['u'][0]
		
		freeEnergies.write("TOTAL: %f \n\n" % (total))
		freeEnergies.close()