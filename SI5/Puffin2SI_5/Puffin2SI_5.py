# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
import numpy as np
import h5py
import tables
import os
import sys
import datetime
now = datetime.datetime.now()
print 'Conversion time: ',now.strftime("%Y-%m-%d %H:%M:%S")
print 'Current directory: ',os.getcwd()

ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=3e+08

if len(sys.argv)==2:
   step_name_in=sys.argv[1]
   print 'Processing file:', step_name_in
else:
   print 'Usage: Puffin2SI_5 <TimeStep Number e.g. 750> \n'
   sys.exit(1)  


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



os.system("sddsxref %s %s %s %s %s %s %s %s"\
%(xdatfile,ydatfile,z2datfile,RE_datfile,IM_datfile,Gamma_datfile,Chi_datfile,outfile))
#os.system("sddscombine -overWrite %s %s %s"\
#%(outfile,ParamDataFile,final_output))

print 'Converting SDDS Puffin output to HDF5...'

final_output_hdf=step_no+'_all.h5'
os.system("sdds2hdf %s %s"\
%(final_output,final_output_hdf))
os.remove(final_output)

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



print 'Processing data...'

Z=Lc_f[0]*Z2_f


X=np.sqrt(Lg_f[0]*Lc_f[0])*X_f[:]
Y=np.sqrt(Lg_f[0]*Lc_f[0])*Y_f[:]

Gamma=Gamma_f[:]*gamma_0_f[0]
P=me*c*((Gamma**2)-1)


Px=RE_PPerp_f[:]*me*c*a_u_f[0]
Py=-1*IM_PPerp_f[:]*me*c*a_u_f[0]
Pz=np.sqrt(P**2-Px**2-Py**2)


NE=s_chi_bar_f*npk_bar[0]

print np.shape(X)
print np.shape(Px)
print np.shape(Y)
print np.shape(Py)
print np.shape(Z)
print np.shape(Pz)
print np.shape(NE)
full_array=np.vstack((X,Px,Y,Py,Z,Pz,NE)).T

FinalHDF_File=step_no+'_PF.si5'

print 'Creating final HDF5 output...'
output_file=tables.open_file(FinalHDF_File,'w')
ParticleGroup=output_file.create_array('/','Particles',full_array)

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
ParticleGroup._v_attrs.FXFELSourceFileOrigin='PUFFIN'
ParticleGroup._v_attrs.FXFELSourceDirectoryName=os.getcwd()

#Close the file
output_file.close()