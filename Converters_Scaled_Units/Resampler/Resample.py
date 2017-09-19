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

Pi=np.pi                    # Pi number taken from 'numpy' as more precise than just 3.1415
k_u=251.327412              # Undulator wave number default=628 k_u=2*Pi/l_w
a_u=1.225699                # undulator parameter ? a_u=a_w
c=3.0e+8                    # Speed of light
m=9.11e-31                  # mass of electron
e_0=8.854E-12               # vacuum permitivity
e_ch=1.602e-19


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
PX  = m*c*scipy.ndimage.interpolation.zoom(Particles[:,1], ResampleRatio,order=0)
PY  = m*c*scipy.ndimage.interpolation.zoom(Particles[:,3], ResampleRatio,order=0)
PZ  = m*c*scipy.ndimage.interpolation.zoom(Particles[:,5], ResampleRatio,order=0)




# Bunch particles to nearest wavelenght (equispace)


size_x=max(X)-min(X)
size_y=max(Y)-min(Y)
size_z=max(Z)-min(Z)
xyz = np.vstack([X,Y,Z]).T
# This is number of bins just to calculate initial data - don't change if not sure
binnumber=100
cube_volume=(size_x*size_y*size_z)/float(binnumber**3)
print 'Volume of sample cube: ',cube_volume
print 'Size of sample X,Y,Z: ', size_x,size_y,size_z
H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=NE.flat)
n_p=float(np.amax(H))/cube_volume
print 'Particle density Np= ',n_p

# Calculate some needed values: Omega,Rho,Lc,Lambda_U,Lambda_R,
p_tot=np.sqrt((PX[:]**2)+(PY[:]**2)+(PZ[:]**2))
gamma=(np.sqrt(1+(p_tot/(m*c))**2))
gamma_0=np.mean(gamma)
omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))
lambda_u=(2*Pi)/k_u

# Calculated for planar undulator -> (1+(a^2)/2)
lambda_r=(lambda_u/(2.0*gamma_0**2.0))*(1+(a_u**2.0)/2.0)
Lc=lambda_r/(4.0*Pi*rho)
print 'Lc = ',Lc
print 'Rho = ',rho
SlicesMultiplyFactor=1
NumberOfSlices=int(SlicesMultiplyFactor*(size_z)/(4.0*Pi*rho*Lc))
print 'Number of slices = ',NumberOfSlices
step=size_z/NumberOfSlices
minz=np.min(Z)
maxz=np.max(Z)
# End of inital data calculations

zgrid=np.arange(minz,maxz,step)
print 'grid shape = ',np.shape(zgrid)
mypoints=Z
kcells=np.searchsorted(zgrid,mypoints[:])-1
zgridded=zgrid[kcells]

# Open output file 
output_file=tables.open_file(file_name_base+'_RES.h5','w')

# Merge all data into one array
CompleteArray = np.vstack([X.flat,(PX/(m*c)).flat,Y.flat,(PY/(m*c)).flat,zgridded.flat,(PZ/(m*c)).flat,NE.flat]).T

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

    