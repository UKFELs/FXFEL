# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
# WARNING !!!!

# This script is draft for converting Puffin data into Elegant
# Elegant needs all particles to be of same charge while Puffin

# The SDDSToolkit is essential !!

import numpy as np
import h5py
import tables
import os
import sys

# Some basic constants
ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=3e+08

# Read the name of the file - as input give the step number you want to use
if len(sys.argv)==2:
   step_name_in=sys.argv[1]
   print 'Processing file:', step_name_in
else:
   print 'Usage: Puf2Ele <TimeStep Number e.g. 750> \n'
   sys.exit(1)  

# Merge Puffin output into one array.
# Find the name of the files
print 'Merging Puffin output files...'
step_no=str(step_name_in)
xdatfile=step_no+'_XDataFile.dat'
ydatfile=step_no+'_YDataFile.dat'
z2datfile=step_no+'_Z2DataFile.dat'
RE_datfile=step_no+'_RE_PPerpDataFile.dat'
IM_datfile=step_no+'_IM_PPerpDataFile.dat'
Gamma_datfile=step_no+'_GammaDataFile.dat'
Chi_datfile='ChiDataFile.dat'
ZDataFile=step_no+'_ZDataFile.dat'
outfile=step_no+'_all.sdds'

final_output=step_no+'_all.sdds'


# Call sddsxref to merge Puffin SDDS formatted output into one large file - takes a while...
# Drink tea, coffee, yerba or whatever you like in case of large files.

os.system("sddsxref %s %s %s %s %s %s %s %s"\
%(xdatfile,ydatfile,z2datfile,RE_datfile,IM_datfile,Gamma_datfile,Chi_datfile,outfile))

# Convert SDDS to HDF5 for easier use in Python
# Python can't access SDDS files directly but can access HDF5 withe ease.

print 'Converting SDDS Puffin output to HDF5...'

final_output_hdf=step_no+'_all.h5'
os.system("sdds2hdf %s %s"\
%(final_output,final_output_hdf))
os.remove(final_output)

# Load data from HDF5 into memory - beware of the size...

print 'Loading data from Puffin HDF5...'

f=h5py.File(final_output_hdf,'r')

X_f=np.array(f['/page1/columns/X'])
Y_f=np.array(f['/page1/columns/Y'])
Z2_f=np.array(f['/page1/columns/Z2'])

Gamma_f=np.array(f['/page1/columns/Gamma'])
IM_PPerp_f=np.array(f['/page1/columns/IM_PPerp'])
RE_PPerp_f=np.array(f['/page1/columns/RE_PPerp'])
s_chi_bar_f=np.array(f['/page1/columns/s_chi_bar'])
f.close()
os.remove(final_output_hdf)

os.system("sdds2hdf ParamDataFile.dat param.h5")

f2=h5py.File('param.h5','r')
Lc_f=np.array(f2['/page1/parameters/Lc'])
Lg_f=np.array(f2['/page1/parameters/Lg'])
a_u_f=np.array(f2['/page1/parameters/aw'])
gamma_0_f=np.array(f2['/page1/parameters/gamma_r'])
npk_bar=np.array(f2['/page1/parameters/npk_bar'])
f2.close()
os.remove('param.h5')

print 'Read constants...:'
print 'Npk_bar= ',npk_bar[0]
print 'Lc = ',Lc_f[0]
print 'Lg = ',Lg_f[0]
print 'aw = ',a_u_f[0]
print 'gamma_r = ',gamma_0_f[0]

os.system("sdds2hdf %s Z.h5"\
%(ZDataFile))
z_file=h5py.File('Z.h5','r')
z_param_f=np.array(z_file['/page1/columns/Z'])
z_file.close()
os.remove('Z.h5')


print 'Z = ',z_param_f[0]

z_bar=z_param_f[0]
print 'z_bar = ', z_bar

total_charge=np.sum(s_chi_bar_f[0])

print 'Total charge= ',total_charge

print 'Processing data...'

# Start to convert coordinates into Elegant format.

z=Lg_f[0]*z_bar

S=(Lc_f[0]*Z2_f[:])+z



X=np.sqrt(Lg_f[0]*Lc_f[0])*X_f[:]
Y=np.sqrt(Lg_f[0]*Lc_f[0])*Y_f[:]

Gamma=Gamma_f[:]*gamma_0_f[0]
P=me*c*((Gamma**2)-1)


Px=RE_PPerp_f[:]*me*c*a_u_f[0]
Py=-1*IM_PPerp_f[:]*me*c*a_u_f[0]
Pz=np.sqrt(P**2-Px**2-Py**2)

Beta_z=Pz/P

xp=Px/Pz
yp=Py/Pz
m_t=z/(c*Beta_z)


# Built array in Elegant format (x,xp,y,yp,t,p)

full_array=np.column_stack((X,xp,Y,yp,m_t,P))

# Open HDF5 output and SDDS output file, HDF5 is used just for plotting
# 
FinalHDF_File=step_no+'_full_array.h5'
FinalSDDS_File=step_no+'_full_array.sdds'

print 'Creating final HDF5 output...'
output_file=tables.open_file(FinalHDF_File,'w')
# Create group in HDF file and remove not necessary metadata.
ParticleGroup=output_file.create_array('/','spatialPositions',full_array)
del ParticleGroup._v_attrs.CLASS
del ParticleGroup._v_attrs.VERSION
del ParticleGroup._v_attrs.FLAVOR
del ParticleGroup._v_attrs.TITLE
output_file.close()
# Convert HDF5 into SDDS
os.system("hdf2sdds %s %s -dataset=spatialPositions"\
%(FinalHDF_File,FinalSDDS_File))

print 'Creating SDDS output...'
# Name columns in SDDS file
os.system("sddsconvert %s -rename=columns,D1=x,D2=xp,D3=y,D4=yp,D5=t,D6=p"\
%(FinalSDDS_File))
# Add total charge which is not true for Elegant... 
os.system("sddsprocess %s -define=parameter,Charge,%s,units=C"\
%(FinalSDDS_File,total_charge))

os.remove(FinalSDDS_File+'~')


# Add metadata to HDF5 file just to use it for plotting in ViSit
print 'Adding some magic...'
output_file=tables.open_file(FinalHDF_File,'a')
ParticleGroup=output_file.root.spatialPositions
output_file=tables.open_file(step_no+'_full_array.h5','a')
boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
boundsGroup._v_attrs.vsType='limits'
boundsGroup._v_attrs.vsKind='Cartesian'
timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'
ParticleGroup._v_attrs.vsType='variableWithMesh'
ParticleGroup._v_attrs.vsTimeGroup='time'
ParticleGroup._v_attrs.vsNumSpatialDims = 3
ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
ParticleGroup._v_attrs.vsLabels='x,xp,y,yp,t,p'
output_file.close()


print 'Created SDDS: ',FinalSDDS_File
print 'Created HDF5: ',FinalHDF_File
#sys.exit(1)