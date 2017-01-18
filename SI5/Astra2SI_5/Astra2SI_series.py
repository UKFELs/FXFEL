# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:55:37 2015

@author: piotrt
"""
import numpy as np
import sys
import tables
import datetime
#set some variables, material constants
now = datetime.datetime.now()
print 'Conversion time: ',now.strftime("%Y-%m-%d %H:%M:%S")

e_ch=1.602e-19

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Astra2SI_5 <FileName> \n'
   sys.exit(1)  
#file_name_base  = (file_name_in.split('.')[0]).strip()
#    ls *.*.001 | xargs -L1 python Astra2SI_5.py
   
fn  = file_name_in.split('.')
print fn[0],fn[1],fn[2]
file_name_base = '{0}_{2}_{0}'.format(fn[0],fn[2],fn[0])
file_name_base = str(fn[0]+fn[2]+'_'+fn[1])
# Open and load file, then read as lines into array
file_in=open(file_name_in,'r')
lines=file_in.readlines()
file_in.close()
data=[]
for line in lines:
		p=line.split()
		data.append(p)
data=[line.split() for line in lines]

# define array as float
data2=np.array(data)
data2=data2.astype(float)
print(data2.shape)

# assign chosen rows from loaded array
x=data2[:,0]
y=data2[:,1]
z=data2[:,2]
px=data2[:,3]
py=data2[:,4]
pz=data2[:,5]
dumb_zero=data2[:,6]
macrocharge=data2[:,7]

# prepare array - needed for Python to work properly
norm_px=px[:]
norm_py=py[:]
norm_pz=pz[:]
norm_x=x[:]
norm_y=y[:]
norm_z=z[:]

# scale first-reference particle
for i in range(1,len(x)):
#    norm_px[i]=(px[i]+px[0])*me
#    norm_py[i]=(py[i]+py[0])*me
#    norm_pz[i]=(pz[i]+pz[0])*me 
#    norm_x[i]=(x[i]+x[0])
#    norm_y[i]=(y[i]+y[0])
#    norm_z[i]=(z[i]+z[0])
    norm_px[i]=(px[i])*5.36E-28
    norm_py[i]=(py[i])*5.36E-28
    norm_pz[i]=(pz[i]+pz[0])*5.36E-28 
    norm_x[i]=(x[i])
    norm_y[i]=(y[i])
    norm_z[i]=(z[i]+z[0])


norm_px[0]=(px[0])*5.36E-28
norm_py[0]=(py[0])*5.36E-28
norm_pz[0]=(pz[0])*5.36E-28 


x_px_y_py_z_pz_NE = np.vstack([norm_x,norm_px,norm_y,norm_py,norm_z,norm_pz,((macrocharge*(-1.E-9))/e_ch)]).T




output_file=tables.open_file(file_name_base+'.h5','w')

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
ParticleGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
ParticleGroup._v_attrs.FXFELSourceFileOrigin='ASTRA'
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
#Close the file
output_file.close()



