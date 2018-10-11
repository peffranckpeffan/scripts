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

def call_subprocess(command, directory, s):
	subp = subprocess.Popen(command, cwd=directory, shell=s)
	subp.wait()

def update_dummy(colvar_dir, colvar_name_list, center_file, selection, stage):
	print("AQUI:"+center_file)
	call_subprocess("env center_file='"+center_file+"' selection='"+selection+"' vmd -dispdev text -e ../scripts/lib-tcl/get-center-mass.tcl", "../"+colvar_dir, True)
	
	centerFile = open("../"+colvar_dir+"/center.tmp")
	center=""
	for line in centerFile:
		center=line
	centerFile.close()
	center=center.split(" ")

	for i in range(0,len(colvar_name_list)):
		lines = open("../"+colvar_dir+"/"+colvar_name_list[i], "r").readlines()
		colv = open("../"+colvar_dir+"/"+colvar_name_list[i], "w")
		update_dummy = 0
		z=center[2].strip()
		for line in lines:
			if ((stage=='equil') and "posit-prot" in line):
				update_dummy = 1

			elif (stage=='smd') and ("posit-XY-lig" in line or "posit-Z-lig" in line):
				update_dummy = 1
				z=str(float(center[2].strip())-4)

			if ("dummyatom" in line) and (update_dummy==1):
				line = "	dummyatom ("+center[0]+", "+center[1]+", "+z+") \n"
				update_dummy=0

			colv.write(line)

		colv.close()

	os.remove("../"+colvar_dir+"/center.tmp")

def createDir(name):
	if not os.path.exists(name):
		os.makedirs(name)
#Copy the lines from a file to another file

def copyFileLines(fileFrom, fileTo):
	file = open(fileFrom)
	for line in file:
		fileTo.write(line)

goFoward=True
stage=""
ifExec=""
try:
	stage=str(sys.argv[1])
except:
	print("ERROR: It's necessary to inform the stage of the analysis.")
	goFoward=False

if goFoward:
	try:
		ifExec=sys.argv[2]
	except:
		print("The stage \""+stage+"\" will be generated, but the simulation will not be perfomed.")
		time.sleep(5)

	if(stage == 'restraints' and ifExec=="exec"):
		try:
			res_folder=str(sys.argv[3])
		except:
			print("In the restraint stage, to run the simulation it's necessary to inform the folder.")
			goFoward=False


if goFoward:
	if (stage == "equil"):

		#Mount the system
		call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-equil.tcl", "../equil", True)

		#Update the dummy atom
		update_dummy('equil', ['colvar.inp'], "../common/solvate.pdb", "segname A and backbone", stage)

		if (ifExec=="exec"):
			file_eq=Path('../equil/out_eq.dcd')
			if (file_eq.exists()):
				call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_eq2 > conf_eq2.log", "../equil", True)
			else:
				call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_eq > conf_eq.log", "../equil", True)
				call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_eq2 > conf_eq2.log", "../equil", True)
			

	if (stage == "smd"):

		#Mount the system
		call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../SMD", True)

		#Update dummyatom
		update_dummy('SMD', ['colv-smd.inp'], "../equil/out_eq2.restart.coor", "segname B and backbone", stage)

		if (ifExec=="exec" and Path('../equil/out_eq.dcd').exists()):
			call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_smd > conf_smd.log", "../SMD", True)

	if (stage == "us"):

		call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../US", True)

		update_dummy('US', ['basic_colv'], "../common/solvate.pdb", "segname A and backbone", 'equil')
		update_dummy('US', ['basic_colv'], "../equil/out_eq2.restart.coor", "segname B and backbone", 'smd')

		#PREPARING MINIMIZATION
		steps=500000
		#initial window in posit_wind + space_bet_wind (3.5)
		posit_wind=3.0
		space_bet_wind=0.5

		for i in range(42):
			if (i < 10):
				window="0"+str(i)
			else:
				window=i
			
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
			copyFileLines("../US/basic_conf", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();

			tcl_file = open("../US/u"+str(window)+"/"+colv_inp, "w")
			copyFileLines("../US/basic_colv", tcl_file)
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
			if (i < 10):
				window="0"+str(i)
			else:
				window=i

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
			copyFileLines("../US/basic_conf", conf_file)
			conf_file.write("run "+ str(steps))
			conf_file.close();

		if ifExec=="exec":
			if(not(Path('../US/out_min.dcd').exists())):
				process=["min", "run"]
			else:
				process=["run"]
			for x in range(len(process)):
				for i in range(0,42):
					if (i < 10):
						window="0"+str(i)
					else:
						window=i
					call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_"+process[x]+"-"+str(window)+" > conf_"+process[x]+"-"+str(window)+".log", "../US/u"+str(window), True)

	#CALCULATE RESTRICTIONS
	if(stage == 'restraints'):

		steps=1000000
		restraints=json.load(open("restraints.json"))["restraints"]
		
		call_subprocess("vmd -dispdev text -e ../scripts/lib-tcl/setup-smd.tcl", "../restraints", True)

		basic_colv_list=['basic_colv-bulk', 'basic_colv-lrmsd', 'basic_colv-orient', 'basic_colv-prmsd', 'basic_colv-trans']
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
				copyFileLines("../restraints/basic_conf", conf_file)
				conf_file.write("run "+ str(steps))
				conf_file.close();


				colv_file = open("../restraints/"+restraints[count_rest]["folder"]+str(nr)+"/"+colv_inp, "w")
				copyFileLines("../restraints/basic_colv-" + restraints[count_rest]["type"], colv_file)
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
						if (i < 10):
							window="0"+str(i)
						else:
							window=i
						call_subprocess("namd2 +p6 +setcpuaffinity +devices 0 +pemap 0-5  conf_rest-"+str(window)+" > conf_rest-"+str(window)+".log", "../restraints/"+folders[x]+str(window), True)






