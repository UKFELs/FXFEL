# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:56:42 2015

@author: piotrt
"""
# Import proper libraries: h5py and numpy, needed to handle HDF5
import h5py
import numpy as np
import sys

m=9.11e-31

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: VSim2MASP <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


# Load file named 36.h5 - other name of file needs change here
f=h5py.File(file_name_in,'r')
# Read the proper table - data from Bernard Hiddings VSim use name 'HeEletrons'
# for array storing X,Y,Z,Px,Py,Pz,0,Weight of particle (that is their format)
# Assign the table 'HeElectrons' to Electrons dataset/table
Electrons=np.array(f['HeElectrons'])
# Assign number of records (macroparticles) as n
n=len(Electrons)
# Create output file named output.txt
out=open(file_name_base+'_2MASP.txt','w')
# Calculate the particle number per record (experimental)
# Assign table 'HeElectrons' to variable dset
dset=f['HeElectrons']
# Read the attribute/variable 'numPtclsInMacro' from HDF5 metadata and assign
# its value to ParticleNumber variable
ParticleNumber=dset.attrs['numPtclsInMacro']
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
for i in range(n):                                  
                out.write("%.8e" %(Electrons[i,0]) + \
                " %.8e" %(Electrons[i,3]*m) + \
                " %.8e" %(Electrons[i,1]) + \
                " %.8e" %(Electrons[i,4]*m) + \
                " %.8e" %(Electrons[i,2]) + \
                " %.8e" %(Electrons[i,5]*m) + \
# The lines below are experimental as it is not decided yet how to handle
# weight of particles from HDF5 source files - a matter of choice which is correct - 3 to choose.
#                " %.8e" %(Electrons[i,7]) + \
# This line below transfers directly weight to ASCII - it works but is from from physical point of view
#                 " %.8e" %(Electrons[i,7]*ParticlePerRecord) + "\n") 
# The line below makes some manipulation of weight of particles (calculates them) but not sure if it is right
# According to Jonny the number of particles should be taken from numPtclsInMacro*weight... Well.    
                 " %.8e" %(Electrons[i,7]*ParticleNumber) + "\n")                
                
#Close the files                
out.close()
f.close()
# End of script