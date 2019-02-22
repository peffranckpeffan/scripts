
### BE SURE THAT THE ENERGY CALCS MATCH THE DISANG FILE ORDER

import numpy as np # numerical array library
import subprocess as sp
import sys,os
from akima import interpolate
import pymbar # multistate Bennett acceptance ratio
import timeseries # timeseries analysis
import math
import pprint as pp

### Arguments
aur = sys.argv[1] # a or u or r or d
temp = float(sys.argv[2]) # temp
try:
  phase = int(sys.argv[3]) # phase
except:
  phase = 0 # phase
kB = 1.381e-23 * 6.022e23 / (4.184 * 1000.0) # Boltzmann constant in kJ/mol/K
beta = 1/(kB * temp) # beta
N_max = 2000000 # Max frames for any simulation window, you should check this if you did some long runs

sys.stdout = open('allsp-'+aur+'.log', 'w')
### Determine Number of umbrellas
K = 0
filename = './'+aur+'%02.0f/restraints.dat' % K
while os.path.isfile(filename):
  K = K+1
  filename = './'+aur+'%02.0f/restraints.dat' % K


R = 1

print  "K= %5.0f  R= %5.0f" % ( K, R ) 


### Allocate storage for simulation data
N = np.zeros([K], np.int32)                       # N_k[k] is the number of snapshots to be used from umbrella simulation k
Neff = np.zeros([K], np.int32)
Nind = np.zeros([K], np.int32)
Nprg = np.zeros([K], np.int32)
rty = ['d']*R                                     # restraint type (distance or angle)
rfc = np.zeros([K,R], np.float64)                 # restraint force constant
rfc2 = np.zeros([K,R], np.float64)                 # restraint force constant
fcmax = np.zeros([R], np.float64)                 # full force constant value used during umbrella portion of work 
req = np.zeros([K,R], np.float64)                 # restraint target value
req2 = np.zeros([K,R], np.float64)                 # restraint target value
val = np.zeros([N_max,K,R], np.float64)           # value of the restrained variable at each frame n
val2 = np.zeros([N_max,K,R], np.float64)           # value of the restrained variable at each frame n
g = np.zeros([K], np.float64)

### Tmp type arrays for energy and spline fitting/integration
u=np.zeros([N_max], np.float64)
x=np.zeros([K],np.float64)
y=np.zeros([K],np.float64)
m=np.zeros([K],np.float64)
s=np.zeros([K],np.float64)

print "Done with array setup\n"


### Read the simulation data
r=0
for k in range(K):
  # Read Equilibrium Value and Force Constant
  if aur == 't':
    with open('./'+aur+'%02.0f/colvar.in' % k, 'r') as f:
      for line in f:
         if 'posit2' in line:                          
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               req[k,r] = float(cols[1])
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc[k,r] = float(cols[1])/2
               break
         if 'posit3' in line:                          
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               req2[k,r] = float(cols[1])
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc2[k,r] = float(cols[1])/2
               break
  elif aur == 'o':
    with open('./'+aur+'%02.0f/colvar.in' % k, 'r') as f:
      for line in f:
         if 'orient2' in line:                          
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               str = cols[1][1:-1]
               req[k,r] = float(str)
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc[k,r] = float(cols[1])/2
               break
  elif aur == 'r' or aur == 'l':
    with open('./'+aur+'%02.0f/colvar.in' % k, 'r') as f:
      for line in f:
         if 'rmsd2' in line:                          
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               req[k,r] = float(cols[1])
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc[k,r] = float(cols[1])/2
               break
  elif aur == 'p' or aur == 'b':
    with open('./'+aur+'%02.0f/colvar.in' % k, 'r') as f:
      for line in f:
         if 'rmsd1' in line:                          
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               req[k,r] = float(cols[1])
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc[k,r] = float(cols[1])/2
               break
  elif aur == 'u':
    with open('./'+aur+'%02.0f/colvar.in' % k, 'r') as f:
      for line in f:
         if 'posit3' in line:
           for line in f:
             cols = line.split()
             if len(cols) != 0 and (cols[0] == "centers"):
               req[k,r] = float(cols[1])
             if len(cols) != 0 and (cols[0] == "forceConstant"):
               rfc[k,r] = float(cols[1])/2
               break
  else:
    sys.exit("not sure about restraint type!")

#req todos os centros dos sistemas nas diferentes janelas
#rfc todos as constantes do sistema nas diferentes janelas
  # pp.pprint(rfc)
  # pp.pprint(req)
  # print('--------------------')
  # Read in Values for restrained variables for each simulation

  filename = './'+aur+'%02.0f/restraints.dat' % k
  infile = open(filename, 'r')
  restdat = infile.readlines()     # slice off first 20 lines  readlines()[20:]
  infile.close()
  # Parse Data
  n = 0
  s = 0
  from_line = -1
  if phase == 0:
    from_line = -1
  for line in restdat:
    s += 1 #so ira analizar o arquivo a partir da linha 500!
    if line[0] != '#' and line[0] != '@' and s > from_line:
      cols = line.split()
      if aur == 'o':
        val[n,k,r] = math.acos(float(cols[2]))
      elif aur == 'u' or aur == 'l':
        val[n,k,r] = float(cols[2])
      elif aur == 't':
         val[n,k,r] = float(cols[1])
         val2[n,k,r] = float(cols[2])
      else:
        val[n,k,r] = float(cols[1])
      n += 1
  N[k] = n
  #val -> valores da respectiva coluna no arquivo restraints.dat
  #sys.exit()
  # Calculate Reduced Potential

  if aur == 'o':
    if rfc[k,0] == 0:
      tmp=np.ones([R],np.float64)*0.001
      u[0:N[k]] = np.sum(beta*tmp[0:R]*((val[0:N[k],k,0:R])**2), axis=1)#->slicing syntax [0:N[k]]
    else:
      u[0:N[k]] = np.sum(beta*rfc[k,0:R]*((val[0:N[k],k,0:R])**2), axis=1)
  elif aur == 't':
    if rfc[k,0] == 0 and rfc2[k,0] != 0:
      tmp=np.ones([R],np.float64)*0.001
      u[0:N[k]] = np.sum(beta*(tmp[0:R]*(((val[0:N[k],k,0:R]-req[k,0:R])**2)) + rfc2[k,0:R]*((val2[0:N[k],k,0:R]-req2[k,0:R])**2)), axis=1) #-> (1/k_bT)*Kx**2
    elif rfc[k,0] != 0 and rfc2[k,0] == 0:
      tmp=np.ones([R],np.float64)*0.001
      u[0:N[k]] = np.sum(beta*(rfc[k,0:R]*(((val[0:N[k],k,0:R]-req[k,0:R])**2)) + tmp[0:R]*((val2[0:N[k],k,0:R]-req2[k,0:R])**2)), axis=1)
    elif rfc[k,0] == 0 and rfc2[k,0] == 0:
      tmp=np.ones([R],np.float64)*0.001
      u[0:N[k]] = np.sum(beta*(tmp[0:R]*(((val[0:N[k],k,0:R]-req[k,0:R])**2)) + tmp[0:R]*((val2[0:N[k],k,0:R]-req2[k,0:R])**2)), axis=1)
    else:
      u[0:N[k]] = np.sum(beta*(rfc[k,0:R]*(((val[0:N[k],k,0:R]-req[k,0:R])**2)) + rfc2[k,0:R]*((val2[0:N[k],k,0:R]-req2[k,0:R])**2)), axis=1)
  else:
    if rfc[k,0] == 0:
      tmp=np.ones([R],np.float64)*0.001
      u[0:N[k]] = np.sum(beta*tmp[0:R]*((val[0:N[k],k,0:R]-req[k,0:R])**2), axis=1) #-> (1/k_bT)*Kx**2
    else:
      u[0:N[k]] = np.sum(beta*rfc[k,0:R]*((val[0:N[k],k,0:R]-req[k,0:R])**2), axis=1)
  
  
#  g[k] = calcg(u[0:N[k]])
#  subs = timeseries.subsampleCorrelatedData(np.zeros([N[k]]),g=g[k])
#  Nind[k] = len(subs)
#  if Nind[k] > 100000:
#    Neff[k] = 100000
#  else:
  Neff[k] = N[k]
  Nind[k] = N[k]


  print  "Processed Window %5.0f.  N= %12.0f.  g= %10.3f   Nind= %12.0f   Neff= %12.0f" % ( k, N[k], g[k], Nind[k], Neff[k] )

print  "Max Neff= %.0f" % ( np.max(Neff) )
Upot = np.zeros([K,K,np.max(Neff)], np.float64)

# Calculate Restraint Energy
for k in range(K):
#  subs = timeseries.subsampleCorrelatedData(np.zeros([N[k]]),g=g[k])
  for l in range(K):
    if aur == 'o':
      Upot[k,l,0:Neff[k]] = np.sum(beta*rfc[l,0:R]*((val[0:Neff[k],k,0:R])**2), axis=1)
    elif aur == 't':
      Upot[k,l,0:Neff[k]] = np.sum(beta*(rfc[l,0:R]*((val[0:Neff[k],k,0:R]-req[l,0:R])**2) + rfc2[l,0:R]*((val2[0:Neff[k],k,0:R]-req2[l,0:R])**2)), axis=1)
    else:
      Upot[k,l,0:Neff[k]] = np.sum(beta*rfc[l,0:R]*((val[0:Neff[k],k,0:R]-req[l,0:R])**2), axis=1)



val=[]

prg = [100]
for p in range(len(prg)):

  Nprg = Neff*prg[p]/100 ## Test integers out only 
  print  "Running MBAR on %.0f percent of the data ... " % ( prg[p] )
  mbar = pymbar.MBAR(Upot, Nprg, verbose = True, method = 'adaptive', initialize = 'BAR')

  print  "Calculate Free Energy Differences Between States"
  [Deltaf, dDeltaf] = mbar.getFreeEnergyDifferences()

  min = np.argmin(Deltaf[0])

  # Write to file
  print  "Free Energy Differences (in units of kcal/mol)"
  print  "%9s %8s %8s %12s %12s" % ('bin', 'f', 'df', 'deq', 'dfc')
  datfile = open('allsp-'+aur+'.%03.0f.dat' % prg[p], 'w')
  for k in range(K):
    if aur == 'r' or aur == 'o' or aur == 'p' or aur == 'b' or aur == 'l':
      print "%10.5f %10.5f %10.5f %12.7f %12.7f" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] )
      datfile.write ( "%10.5f %10.5f %10.5f %12.7f %12.7f\n" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ) )
    elif aur == 't':
      print "%10.5f %10.5f %10.5f %12.7f %12.7f %12.7f %12.7f" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], req2[k,0], rfc[k,0] ,rfc2[k,0])
      datfile.write ( "%10.5f %10.5f %10.5f %12.7f %12.7f\n" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ) )
    elif aur == 'd':
      print "%9.0f %10.5f %10.5f %12.7f %12.7f" % ( k, Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0]/rfc[-1,0] )
      datfile.write ( "%9.0f %10.5f %10.5f %12.7f %12.7f\n" % ( k, Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0]/rfc[-1,0] ) )
    else: # 'u'
      print "%10.5f %10.5f %10.5f %12.7f %12.7f" % ( req[k,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] )
      datfile.write ( "%10.5f %10.5f %10.5f %12.7f %12.7f\n" % ( req[k,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ) )
  datfile.close()
  print "\n\n"

