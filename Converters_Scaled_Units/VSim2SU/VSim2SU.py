# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:56:42 2015

@author: piotrt
"""
# Import proper libraries: h5py and numpy, needed to handle HDF5
#import h5py
import tables
import numpy as np
import sys

c=3.0e+8
m=9.11e-31

if len(sys.argv)==3:
   file_name_in=sys.argv[1]
   array_name_in=sys.argv[2]
   print 'Processing file:', file_name_in,' Array name is: ',array_name_in
else:
   print 'Usage: VSim2SI_5 <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()



f=tables.open_file(file_name_in,'r')

Electrons=f.root._v_children[array_name_in].read()
# Assign number of records (macroparticles) as n
n=len(Electrons)

# Read the attribute/variable 'numPtclsInMacro' from HDF5 metadata and assign
# its value to ParticleNumber variable


ParticleNumber=f.root._v_children[array_name_in].attrs.numPtclsInMacro
# Assign RecordsNumber as n (number of macroparticles in source file)

RecordsNumber=n
# Simple calculate the number of Particles per Record - NOT sure if it is all right
ParticlePerRecord=ParticleNumber/RecordsNumber


# The Axis Z and X are switched - VSim assumes particles travel along X-axis
# while in Puffin we need particles to travel along Z-axis
# X<->Z
# Px<->Pz
# Assign data to array and convert to scaled units (p = p/mc)

x=Electrons[:,2]
px=Electrons[:,5]/c
y=Electrons[:,1]
py=Electrons[:,4]/c
z=Electrons[:,0]
pz=Electrons[:,3]/c
NE=Electrons[:,7]*ParticleNumber

x_px_y_py_z_pz_NE = np.vstack((x,px,y,py,z,pz,NE)).T




output_file=tables.open_file(file_name_base+'_VS.h5','w')

# Create hdf5 file


# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','Particles',x_px_y_py_z_pz_NE)

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
#Close the file
output_file.close()