# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:03:14 2015

@author: piotrt
"""
# The script converts MASP 3D output into 3D Puffin input.
# The MASP format is:
# x,px,y,py,z,pz, weight_of_particle

# Import necessary libraries
import tables
import time
import numpy as np
import csv
import sys

# Reade the name of the input file
if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: MASP2Puf <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


# time - just for benchmarking, not necessary
start = time.time()
def elapsed():
    return time.time() - start

# count data rows, to preallocate array
# Change the file name - twice ! (the next is repeated in line 37, this is for debug reason)
linecount = sum(1 for line in open(file_name_in))

print '\n%.3fs: File has %s rows' % (elapsed(), linecount)

# pre-allocate array and load data into array

m_X = np.zeros((linecount,1), dtype=np.float64)
m_PX = np.zeros((linecount,1), dtype=np.float64)
m_Y = np.zeros((linecount,1), dtype=np.float64)
m_PY = np.zeros((linecount,1), dtype=np.float64)
m_Z = np.zeros((linecount,1), dtype=np.float64)
m_PZ = np.zeros((linecount,1), dtype=np.float64)
m_WGHT = np.zeros((linecount,1), dtype=np.float64)

# Read the file line by line but each row goes to separate array (so you need to know structure it is not good)
# The delimiter is important as MASP produces output with comma while other programs can use space - check before run
# Check the filename - has to be same as the one used previously

f = csv.reader(open(file_name_in, 'r'),delimiter=' ')
for i, row in enumerate(f):
     m_X[i] = float(row[0]),
     m_PX[i] = float(row[1]),
     m_Y[i] = float(row[2]),
     m_PY[i] = float(row[3]),
     m_Z[i] = float(row[4]),
     m_PZ[i] = float(row[5]),
     m_WGHT[i] = float(row[6]),

# Just the time for benchmark
print '%.3fs: loaded' % elapsed()

# Load trump file proper for input file you use
# It is needed to find maximum value of Np

print 'Loading trump.txt...'
tcol1,tcol2,tcol3,n_p_column=np.loadtxt('trump.txt',unpack=True)
n_p=np.max(n_p_column)
print 'Np= ',n_p

# Puffin variables - may need change
e_ch=1.602e-19
c=3.0e+8
m=9.11e-31
Pi=3.1415
k_u=628                   # Undulator wave number
a_u=1                     # undulator parameter ? a_u=a_w
e_0=8.854*10**(-12)       # vacuum permitivity

# Start to calculate data needed for Puffin

p_tot=np.sqrt((m_PX[:]**2)+(m_PY[:]**2)+(m_PZ[:]**2))
gamma=(np.sqrt(1+(p_tot/(m*c))**2))

omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
gamma_0=np.mean(gamma)
print 'Gamma= ',gamma_0

rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2/3))
lambda_u=(2*Pi)/k_u
lambda_r=(lambda_u/(2*gamma_0**2))*(1+a_u**2)

print 'Rho= ', rho
print 'Lambda_u= ',lambda_u
print 'Lambda_r= ',lambda_r

Lc=lambda_r/(4*Pi*rho)
Lg=lambda_u/(4*Pi*rho)

# Puffin output arrays:
z2=m_Z/Lc

x_bar=m_X[:]/(np.sqrt(Lg*Lc))
y_bar=m_Y[:]/(np.sqrt(Lg*Lc))
px_bar=m_PX[:]/(m*c*a_u)
py_bar=m_PY[:]/(m*c*a_u)
Ne=m_WGHT[:]/n_p          # weight


# Combine all read arrays into one
m=np.hstack((x_bar,y_bar,px_bar,py_bar,gamma,z2,Ne))

output_file=tables.open_file(file_name_base+'_Puffin.hdf','w')

# Create hdf5 file


# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','spatialPositions', m)

#Create metadata - currently for Visit to make scatter plots
boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
boundsGroup._v_attrs.vsType='limits'
boundsGroup._v_attrs.vsKind='Cartesian'
timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'
ParticleGroup._v_attrs.vsType='variableWithMesh'
ParticleGroup._v_attrs.vsTimeGroup='time'
ParticleGroup._v_attrs.vsNumSpatialDims = 3
ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
ParticleGroup._v_attrs.vsLabels='x_bar,y_bar,px_bar,py_bar,gamma,z2,Ne'

#Close the file
output_file.close()
# Create Puffin ASCII input file
out_txt=open(file_name_base+'_Puffin.txt','w')

n=len(m_X)

out_txt.write('__[x_bar]_______[y_bar]_______[px_bar]______[py_bar]_____ [Gamma]________[Z2]__________[Ne]_____'+'\n')
for i in range(n): 
        out_txt.write(      "%.6e" %(x_bar[i]) + \
                            " % .6e" %(y_bar[i]) + \
                            " % .6e" %(px_bar[i]) + \
                            " % .6e" %(py_bar[i]) + \
                            " % .6e" %(gamma[i]) + \
                            " % .6e" %(z2[i]) + \
                            " % .6e" %(Ne[i]) + "\n")
out_txt.close()








