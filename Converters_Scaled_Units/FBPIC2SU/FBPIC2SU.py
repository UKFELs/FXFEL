# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:02:41 2017

@author: ptracz
"""
import tables
import h5py
import numpy as np
import sys

m=9.11e-31
c=2.99792458E8

if len(sys.argv)==3:
   file_name_in=sys.argv[1]
   species=sys.argv[2]
   print 'Processing file:', file_name_in
   f=h5py.File(file_name_in,'r')
elif len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
   f=h5py.File(file_name_in,'r')
   root_dir=str(f.keys()[0])
   time_step=str(f[root_dir].keys()[0])
   species_avail = f[root_dir+'/'+time_step+'/particles/'].keys()
   print 'Usage: python FBIC2SU.py <FileName> <ParticleSpecies>\n'
   print 'Available species are:'
   for i in species_avail:
      print i.encode('utf-8')+'\n'
   sys.exit(1)
else:
   print 'Usage: python FBIC2SU.py <FileName> <ParticleSpecies>\n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


root_dir=str(f.keys()[0])
time_step=str(f[root_dir].keys()[0])

m_X=np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/position/x'])
m_Y=np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/position/y'])
m_Z=np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/position/z'])
m_PX=(m*c)*np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/momentum/x'])
m_PY=(m*c)*np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/momentum/y'])
m_PZ=(m*c)*np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/momentum/z'])
m_NE=np.array(f[root_dir+'/'+time_step+'/particles/'+species+'/weighting'])


m_Arr=np.vstack((m_X,m_PX,m_Y,m_PY,m_Z,m_PZ,m_NE)).T

output_file=tables.open_file(file_name_base+'_FBIC.h5','w')


ParticleGroup=output_file.create_array('/','Particles', m_Arr)
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
output_file.close()
