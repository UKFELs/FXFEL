# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
import numpy as np
import h5py
import tables
import os
step_no=str(750)
xdatfile=step_no+'_XDataFile.dat'
ydatfile=step_no+'_YDataFile.dat'
z2datfile=step_no+'_Z2DataFile.dat'
RE_datfile=step_no+'_RE_PPerpDataFile.dat'
IM_datfile=step_no+'_IM_PPerpDataFile.dat'
Gamma_datfile=step_no+'_GammaDataFile.dat'
Chi_datfile='ChiDataFile.dat'
outfile='all.sdds'

ParamDataFile='ParamDataFile.dat'
final_output='all2.sdds'


os.system("sddsxref %s %s %s %s %s %s %s %s"\
%(xdatfile,ydatfile,z2datfile,RE_datfile,IM_datfile,Gamma_datfile,Chi_datfile,outfile))
os.system("sddscombine %s %s %s"\
%(outfile,ParamDataFile,final_output))

final_output_hdf='all2.h5'
os.system("sdds2hdf %s %s"\
%(final_output,final_output_hdf))


f=h5py.File('all2.h5','r')

X_f=np.array(f['/page1/columns/X'])
Y_f=np.array(f['/page1/columns/Y'])
Z2_f=np.array(f['/page1/columns/Z2'])

Gamma_f=np.array(f['/page1/columns/Gamma'])
IM_PPerp_f=np.array(f['/page1/columns/IM_PPerp'])
RE_PPerp_f=np.array(f['/page1/columns/RE_PPerp'])
s_chi_bar_f=np.array(f['/page1/columns/s_chi_bar'])

me=9.11E-31
c=3E+08

#params_file=h5py.File('param.h5','r')
Lc_f=np.array(f['/page2/parameters/Lc'])
Lg_f=np.array(f['/page2/parameters/Lg'])
a_u_f=np.array(f['/page2/parameters/aw'])
gamma_0_f=np.array(f['/page2/parameters/gamma_r'])
print 'Lc = ',Lc_f[0]
print 'Lg = ',Lg_f[0]
print 'aw = ',a_u_f[0]
print 'gamma_r = ',gamma_0_f[0]
z_file=h5py.File('Z.h5','r')
z_param_f=np.array(z_file['/page1/columns/Z'])
t=z_param_f[0]/c
print 'Z = ',z_param_f[0]

z_bar=z_param_f[0]
print 'z_bar = ', z_bar
z=Lg_f[0]*z_bar

S=(Lc_f[0]*Z2_f[:])+z



X=np.sqrt(Lg_f[0]*Lc_f[0])*X_f[:]
Y=np.sqrt(Lg_f[0]*Lc_f[0])*Y_f[:]


Px=RE_PPerp_f[:]*me*c*a_u_f[0]
Py=-1*IM_PPerp_f[:]*me*c*a_u_f[0]

Gamma=Gamma_f[:]*gamma_0_f[0]

full_array=np.column_stack((X,Y,S,Px,Py,Gamma,s_chi_bar_f))
output_file=tables.open_file(step_no+'_full_array.h5','w')

ParticleGroup=output_file.create_array('/','spatialPositions',full_array)
boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
boundsGroup._v_attrs.vsType='limits'
boundsGroup._v_attrs.vsKind='Cartesian'
timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'
ParticleGroup._v_attrs.vsType='variableWithMesh'
ParticleGroup._v_attrs.vsTimeGroup='time'
ParticleGroup._v_attrs.vsNumSpatialDims = 3
ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
ParticleGroup._v_attrs.vsLabels='x,y,s,Px,Py,Gamma,Chi'
output_file.close()

os.system("hdf2sdds 750_full_array.h5 full_array.sdds -dataset=spatialPositions")
os.system("sddsconvert full_array.sdds -rename=columns,D1=x,D2=y,D3=s,D4=Px,D5=Py,D6=Gamma,D7=Chi")
os.remove("full_array.sdds~")
