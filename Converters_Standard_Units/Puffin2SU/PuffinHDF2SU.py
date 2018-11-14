# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:03:14 2015

@author: piotrt
"""
# Import necessary libraries
import tables

import numpy as np
import sys


if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Puffin2SU <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


f=tables.open_file(file_name_in,'r')
print f.root.runInfo._v_attrs.Lg
Lg=f.root.runInfo._v_attrs.Lg
Lc=f.root.runInfo._v_attrs.Lc
a_u=f.root.runInfo._v_attrs.aw
gamma_Scale=f.root.runInfo._v_attrs.gamma_r
npk=f.root.runInfo._v_attrs.npk_bar



Electrons=f.root.electrons.read()

c=3.0e+8
m=9.11e-31


# Assign the data to array and descale the units from p/mc to SI !

m_X = Electrons[:,0]
m_Y = Electrons[:,1]
m_Z = Electrons[:,2]
m_PX = Electrons[:,3]
m_PY = Electrons[:,4]
m_GAMMA = Electrons[:,5]
m_WGHT = Electrons[:,6]






x_su=m_X*np.sqrt(Lg*Lc)
y_su=m_Y*np.sqrt(Lg*Lc)
z_su=m_Z*Lc
px_su=m*c*a_u*m_PX
py_su=-1.0*m*c*a_u*m_PY
gamma=m_GAMMA*gamma_Scale
pz_su=np.sqrt(((m*c)**2.0)*(gamma**2.0-1.0)-px_su**2.0-py_su**2.0)

WGHT=m_WGHT*(npk)
#print np.size(p_tot),p_tot1
FullArray = np.vstack((x_su,px_su/(m*c),y_su,py_su/(m*c),z_su,pz_su/(m*c),WGHT)).T




output_file=tables.open_file(file_name_base+'_SU.h5','w')

# Create hdf5 file


# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','Particles',FullArray)

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







