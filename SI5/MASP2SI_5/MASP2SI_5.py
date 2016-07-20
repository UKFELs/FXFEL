# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:55:37 2015

@author: piotrt
"""
import numpy as np
import tables
import sys
import csv
import time
import datetime
now = datetime.datetime.now()
print 'Conversion time: ',now.strftime("%Y-%m-%d %H:%M:%S")

#set some variables, material constants

c=3.0e+8
e_ch=1.602e-19

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: MASP2SI_5 <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()



# Open and load file, then read as lines into array

start = time.time()
def elapsed():
    return time.time() - start

# count data rows, to preallocate array
# Change the file name - twice ! (the next is repeated in line 37, this is for debug reason)
linecount = sum(1 for line in open(file_name_in))

print '\n%.3fs: File has %s rows' % (elapsed(), linecount)

# pre-allocate array and load data into array

x = np.zeros((linecount,1), dtype=np.float64)
px = np.zeros((linecount,1), dtype=np.float64)
y = np.zeros((linecount,1), dtype=np.float64)
py = np.zeros((linecount,1), dtype=np.float64)
z = np.zeros((linecount,1), dtype=np.float64)
pz = np.zeros((linecount,1), dtype=np.float64)
NE = np.zeros((linecount,1), dtype=np.float64)

# Read the file line by line but each row goes to separate array (so you need to know structure it is not good)
# The delimiter is important as MASP produces output with comma while other programs can use space - check before run
# Check the filename - has to be same as the one used previously

f = csv.reader(open(file_name_in, 'r'),delimiter=',')
for i, row in enumerate(f):
     x[i] = float(row[0]),
     px[i] = float(row[1]),
     y[i] = float(row[2]),
     py[i] = float(row[3]),
     z[i] = float(row[4]),
     pz[i] = float(row[5]),
     NE[i] = float(row[6]),

# Just the time for benchmark
print '%.3fs: loaded' % elapsed()


x_px_y_py_z_pz_NE = np.vstack([x,px,y,py,z,pz,NE]).T




output_file=tables.open_file(file_name_base+'.si5','w')

# Create hdf5 file


# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','Electrons',x_px_y_py_z_pz_NE)

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
ParticleGroup._v_attrs.vsLabels='x,px,y,py,z,pz,NE'
ParticleGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
ParticleGroup._v_attrs.FXFELSourceFileOrigin='MASP'
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
#Close the file
output_file.close()



