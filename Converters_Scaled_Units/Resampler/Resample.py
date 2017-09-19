# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 13:04:19 2017

@author: ptracz
"""

import numpy as np
import tables
import sys
import time
import datetime
import scipy.ndimage.interpolation

now = datetime.datetime.now()
print 'Current time: ',now.strftime("%Y-%m-%d %H:%M:%S")



if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Resample <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()

 
# time - just for benchmarking,
start = time.time()
def elapsed():
    return time.time() - start


FileIn=tables.open_file(file_name_in,'r')
Particles=FileIn.root.Particles.read()
Particles=Particles[Particles[:,4].argsort()]
#*************************************************************
# Density factor i.e multiplier for number of particles
DesiredOutputParticles=10.E6
# Fix the output number of macroparticles to DesiredOutputParticles
DensityFactor=DesiredOutputParticles/len(Particles)
print 'Density Factor = ',DensityFactor
#*************************************************************

ResampleRatio=DensityFactor

X  = scipy.ndimage.interpolation.zoom(Particles[:,0], ResampleRatio,order=0)
Y  = scipy.ndimage.interpolation.zoom(Particles[:,2], ResampleRatio,order=0)
Z  = scipy.ndimage.interpolation.zoom(Particles[:,4], ResampleRatio,order=0)
NE = (scipy.ndimage.interpolation.zoom(Particles[:,6], ResampleRatio,order=0))/ResampleRatio
PX  = scipy.ndimage.interpolation.zoom(Particles[:,1], ResampleRatio,order=0)
PY  = scipy.ndimage.interpolation.zoom(Particles[:,3], ResampleRatio,order=0)
PZ  = scipy.ndimage.interpolation.zoom(Particles[:,5], ResampleRatio,order=0)



# Open output file 
output_file=tables.open_file(file_name_base+'_RES.h5','w')

# Merge all data into one array
CompleteArray = np.vstack([X.flat,PX.flat,Y.flat,PY.flat,Z.flat,PZ.flat,NE.flat]).T

print 'Saving the output to files...'

# Create Group in HDF5 and add metadata for Visit
ParticleGroup=output_file.create_array('/','Particles',CompleteArray)
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
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
ParticleGroup._v_attrs.FXFELDensityFactor=DensityFactor
#Close the file
output_file.close()

    