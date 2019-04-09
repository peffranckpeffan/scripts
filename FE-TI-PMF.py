
### BE SURE THAT THE ENERGY CALCS MATCH THE DISANG FILE ORDER

import math
import numpy as np # numerical array library
import sys
import os
import timeseries # timeseries analysis
from akima import interpolate

### Arguments
aur = sys.argv[1] # a or u or r
N_max = 2000000    

### Determine Number of umbrellas
K = 0
filename = './'+aur+'%02.0f/restraints.dat' % K
while os.path.isfile(filename):
  K = K+1
  filename = './'+aur+'%02.0f/restraints.dat' % K

R=1

### Functions
def factors(n):   # Return list of integer factors
  facs=[]
  sqrt=int(round(np.sqrt(n)+0.5))
  i = 1
  while i <= sqrt:
    if n % i == 0:
      facs.append(i)
      j=n/i
      if j != i:
        facs.append(j)
    i += 1
  return sorted(facs, key=int)

def nearestmax(n):  # Nearest number that is 100 less than input and divisible by two
  int=[]
  numfac=[]
  maxfac=0
  if n % 2 == 0:
    low=n-100
    high=n
  else:
    low=n-101
    high=n-1
  if low < 0:
    low=0
  for i in range(low,high+2,2):
    numfac = len(factors(i))
    if numfac >= maxfac:
      maxfac = numfac
      mostfac = i
  return mostfac

def seom(N,arr):
  Facs=[]
  SEOM=[]
  Facs=factors(N)
  Bn=np.zeros([len(Facs)], np.int32) # Number of blocks for a given block size
  Bmean=np.zeros([len(Facs),Facs[-1]], np.float64) # Means for each block of each block size
  SEOM=np.zeros([len(Facs)-2],np.float64)
  for i in range(len(Facs)-2):   # Run over all block sizes except the final two: two blocks, one block.
    for j in range(Facs[-i-1]):  # Run over all blocks in the data for a specific size.
      Bmean[i,j]=np.mean(arr[j*Facs[i]:(j+1)*Facs[i]])
    Bn[i]=j+1
    SEOM[i]=np.std(Bmean[i,0:Bn[i]],ddof=0)/np.sqrt(Bn[i]-1)
  return np.max(SEOM)
  #hist,edges = np.histogram(SEOM,range=[np.min(SEOM),np.max(SEOM)],bins=50,)  ### edges = len(hist)+1
  #return np.max(SEOM),edges[np.argmax(hist)+1]  # Assume the blocking plateau corresponds to the max histogram count!

### Allocate storage for simulation data
N = np.zeros([K], np.int32)                       # N_k[k] is the number of snapshots from umbrella simulation k
u = np.zeros([K,N_max], np.float64)               # Forces, as in translation = dChemPotential/dr or attachment = dChemPotential/dlamda 
sti = np.zeros([K], np.float64)                   # Statistical Inefficiency
rty = ['d']*R                                     # restraint type (distance or angle)
rfc = np.zeros([K,R], np.float64)                 # restraint force constant
fcmax = np.zeros([R], np.float64)                 # full force constant value used during umbrella portion of work 
req = np.zeros([K,R], np.float64)                 # restraint target value
val = np.zeros([N_max,K,R], np.float64)           # value of the restrained variable at each frame n

### Tmp type arrays for spline fitting/integration
x=np.zeros([K],np.float64)
y=np.zeros([K],np.float64)
m=np.zeros([K],np.float64)
s=np.zeros([K],np.float64)

### Read the simulation data
r=0
for k in range(K):
  # Read Equilibrium Value and Force Constant
  if aur == 'a':
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
  elif aur == 'r':
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
               #print rfc[k,r]
               break
  else:
    sys.exit("not sure about restraint type!")

  # Read in Values for restrained variables for each simulation
  filename = './'+aur+'%02.0f/restraints.dat' % k
  infile = open(filename, 'r')
  restdat = infile.readlines()     # slice off first 20 lines  readlines()[20:]
  infile.close()
  # Parse Data
  n = 0
  t = 0
  for line in restdat:
    t += 1
    if line[0] != '#' and line[0] != '@' and t > 200:
      cols = line.split()
      if aur == 'o':
        val[n,k,r] = math.acos(float(cols[2]))
      elif aur == 'u':
        val[n,k,r] = float(cols[2])
      else:
        val[n,k,r] = float(cols[1])
      n += 1
  N[k] = n


### Identify the Max Force Constant Value (ie, the "attached, bound state")
for r in range(R):
  fcmax[r] = np.max(rfc[0:K,r])

print(N)
prg = [100]
for p in range(len(prg)):
  # Calculate Restraint Energy and Write Forces
  datfile = open('forces-'+aur+'.%03.0f.dat' % prg[p], 'w')
  for k in range(K):
    Nprg = nearestmax(N[k]*prg[p]/100)

    if aur == 'o': ### Attach Orientational
      x[k] = rfc[k,0]/fcmax[0]
      u[k,0:Nprg] = np.sum(fcmax[0:R]*((val[0:Nprg,k,0:R])**2), axis=1)
    elif aur == 'u': ### Umbrella
      x[k] = req[k,0]
      u[k,0:Nprg] = (2.0*rfc[k,0]*(val[0:Nprg,k,0]-req[k,0]))
    else: ### Attach Restraints
      x[k] = rfc[k,0]/fcmax[0]
      u[k,0:Nprg] = np.sum(fcmax[0:R]*((val[0:Nprg,k,0:R]-req[k])**2), axis=1)

 
    m[k] = np.mean(u[k,0:Nprg])
    s[k] = seom(Nprg,u[k,0:Nprg]) ## Blocking Method
    datfile.write( "%15.7f %15.7f %15.7f\n" % (x[k], m[k], s[k])) #, Nprg, np.std(u[k,0:Nprg])
  datfile.close()
  
  ### Prepare to integrate: Make an array to hold spline (xspl,yspl) and record indices for simulation points (idx)
  BootCyc = 10000  # Num of Boot strap cycles
  intg = np.zeros( [K, BootCyc], np.float64 ) # Integration value at each k for each boot cycle
  xspl = np.zeros( [0], np.float64 ) # array for x dimension spline points
  idx = np.zeros( [K], np.int32 ) # index location of the k window in the spline array
  idx[0] = 0
  for k in range(1,K):
    xspl = np.append(xspl,np.linspace(x[k-1],x[k],num=100,endpoint=False)) # 100 spline poins between windows
    idx[k]=len(xspl)
  xspl = np.append(xspl,x[-1])
  
  ### Integrate with bootstrapping
  for i in range(BootCyc):
    for k in range(K):
      y[k]=np.random.normal(m[k],s[k],1)
    yspl = interpolate(x,y,xspl)
    for k in range(K):
      if aur == 'a' or aur == 'o' or aur == 'r': ### Attach/Release
        intg[k,i] = np.trapz(yspl[0:idx[k]],xspl[0:idx[k]]) # not a negative because it's work by the system, not external.  test with w = -fd
      else: # Umbrella Translate
        intg[k,i]=-1*np.trapz(yspl[0:idx[k]],xspl[0:idx[k]])
  
  ### Write Integration
  datfile = open('int-'+aur+'.%03.0f.dat' % prg[p], 'w')
  for k in range(K):
    datfile.write("%15.7f %15.7f %15.7f\n" % ( x[k], np.mean(intg[k]), np.std(intg[k]) ) )
  datfile.close()
  
  ### Write Spline
  yspl = interpolate(x,m,xspl)
  datfile = open('spline-'+aur+'.%03.0f.dat' % prg[p], 'w')
  for i in range(yspl.size):
    datfile.write("%15.7f %15.7f\n" % (xspl[i],yspl[i]) )
  datfile.close()
  
