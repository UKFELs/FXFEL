# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:55:37 2015

@author: piotrt
"""
import numpy as np
import tables
import sys
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


data=[]
with open(file_name_in, 'r') as file_in:
    for line in file_in:
            p=line.split()
            data.append(p)
            data=[line.split() for line in file_in]


# define array as float
data2=np.array(data,dtype=np.float64)
#data2=data2.astype(float)
print(data2.shape)

# assign chosen rows from loaded array
x=data2[:,0]
y=data2[:,2]
z=data2[:,4]
px=data2[:,1]
py=data2[:,3]
pz=data2[:,5]
NE=data2[:,6]



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



