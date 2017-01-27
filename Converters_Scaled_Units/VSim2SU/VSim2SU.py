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

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: VSim2SI_5 <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


# Load file named 36.h5 - other name of file needs change here
#f=h5py.File(file_name_in,'r')
f=tables.open_file(file_name_in,'r')
# Read the proper table - data from Bernard Hiddings VSim use name 'HeEletrons'
# for array storing X,Y,Z,Px,Py,Pz,0,Weight of particle (that is their format)
# Assign the table 'HeElectrons' to Electrons dataset/table
Electrons=f.root.HeElectrons.read()
# Assign number of records (macroparticles) as n
n=len(Electrons)
# Create output file named output.txt
out=open(file_name_base+'.si5','w')
# Calculate the particle number per record (experimental)
# Assign table 'HeElectrons' to variable dset
#dset=f['HeElectrons']
# Read the attribute/variable 'numPtclsInMacro' from HDF5 metadata and assign
# its value to ParticleNumber variable
#ParticleNumber=dset.attrs['numPtclsInMacro']

ParticleNumber=f.root.HeElectrons.attrs.numPtclsInMacro
# Assign RecordsNumber as n (number of macroparticles in source file)
RecordsNumber=n
# Simple calculate the number of Particles per Record - NOT sure if it is all right
ParticlePerRecord=ParticleNumber/RecordsNumber
# Loop has defualt format (what is in - goes out) NOT good for MASP - therefore commented and ignored by compiler
#for i in range(n):                                  
#                out.write("%.8e" %(Electrons[i,0]) + \
#                "\t%.8e" %(Electrons[i,1]) + \
#                "\t%.8e" %(Electrons[i,2]) + \
#                "\t%.8e" %(Electrons[i,3]) + \
#                "\t%.8e" %(Electrons[i,4]) + \
#                "\t%.8e" %(Electrons[i,5]) + \
#                "\t%.8e" %(Electrons[i,6]) + \
#                "\t%.8e" %(Electrons[i,7]) + "\n")

# Loop adapted to MASP file format: X,Px,Y,Py,Z,Pz, 7 is "something"


# Write the whole array in MASP input format X,Px,Y,Py,Z,Pz,Weight_of_Particle
# This is done in loop - reading each item from Electrons array
# !!!! WARNING !!!!
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