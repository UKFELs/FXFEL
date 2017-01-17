# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
# Experimental code - visualises EM fields from Puffin data

import numpy as np
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
   print 'Usage: Fields2SI_5 <TimeStep Number e.g. 750> \n'
   sys.exit(1)  
step_no=str(step_name_in)


IM_SDDS=step_no+'_IM_ADataFile.dat'
RE_SDDS=step_no+'_RE_ADataFile.dat'

final_output_hdf_IM=step_no+'_IM_ADataFile.h5'
final_output_hdf_RE=step_no+'_RE_ADataFile.h5'

os.system("sdds2hdf %s %s"\
%(IM_SDDS,final_output_hdf_IM))
os.system("sdds2hdf %s %s"\
%(RE_SDDS,final_output_hdf_RE))







print 'Merging Puffin output files...'
step_no=str(step_name_in)
xdatfile=step_no+'_XDataFile.dat'
ydatfile=step_no+'_YDataFile.dat'
z2datfile=step_no+'_Z2DataFile.dat'
outfile=step_no+'_all.sdds'

os.system("sddsxref %s %s %s %s"\
%(xdatfile,ydatfile,z2datfile,outfile))
#os.system("sddscombine -overWrite %s %s %s"\
#%(outfile,ParamDataFile,final_output))

print 'Converting SDDS Puffin output to HDF5...'

final_output_hdf=step_no+'_POS.h5'
os.system("sdds2hdf %s %s"\
%(outfile,final_output_hdf))



f=tables.open_file(final_output_hdf,'r')
X_f=f.root.page1.columns.X.read()
Y_f=f.root.page1.columns.Y.read()
Z2_f=f.root.page1.columns.Z2.read()



print 'Loading data from Puffin HDF5...'

#f=h5py.File(final_output_hdf,'r')
f_IM=tables.open_file(final_output_hdf_IM,'r')
f_RE=tables.open_file(final_output_hdf_RE,'r')
IM_A=f_IM.root.page1.columns.IM_A.read()
RE_A=f_RE.root.page1.columns.RE_A.read()

os.system("sdds2hdf ParamDataFile.dat ParamDataFile.h5")

f2=tables.open_file('ParamDataFile.h5','r')
nX=f2.root.page1.parameters.nX.read()
nY=f2.root.page1.parameters.nY.read()
nZ2=f2.root.page1.parameters.nZ2.read()

Lc_f=f2.root.page1.parameters.Lc.read()
Lg_f=f2.root.page1.parameters.Lg.read()

print 'Read constants...:'
print 'nX = ',nX[0]
print 'nY = ',nY[0]
print 'nZ2 = ',nZ2[0]

Z=Lc_f[0]*Z2_f


X=np.sqrt(Lg_f[0]*Lc_f[0])*X_f[:]
Y=np.sqrt(Lg_f[0]*Lc_f[0])*Y_f[:]



IM_A_RSHP=np.reshape(IM_A,(nX[0],nY[0],nZ2[0]),order='F')
RE_A_RSHP=np.reshape(RE_A,(nX[0],nY[0],nZ2[0]),order='F')



print 'Processing data...'


#full_array=np.vstack((X,Px,Y,Py,Z,Pz,NE)).T

FinalHDF_File=step_no+'_Fields.h5'

print 'Creating final HDF5 output...'
output_file=tables.open_file(FinalHDF_File,'w')

FULL_FIELD=IM_A_RSHP**2+RE_A_RSHP**2

FieldGroup=output_file.create_array('/','Fields',FULL_FIELD)

#Create metadata 

boundsGroup=output_file.create_group('/','globalLimits','')
boundsGroup._v_attrs.vsKind='Cartesian'

mins=np.array([min(X),min(Y),min(Z)])
maxs=np.array([max(X),max(Y),max(Z)])

#boundsGroup._v_attrs.vsLowerBounds=np.array([np.min(FULL_FIELD[1,:,:]),np.min(FULL_FIELD[:,1,:]),np.min(FULL_FIELD[:,:,1])])
boundsGroup._v_attrs.vsLowerBounds=mins
boundsGroup._v_attrs.vsType='limits'
#boundsGroup._v_attrs.vsUpperBounds=np.array([np.max(FULL_FIELD[1,:,:]),np.max(FULL_FIELD[:,1,:]),np.max(FULL_FIELD[:,:,1])])
boundsGroup._v_attrs.vsUpperBounds=maxs
boundsGroup._v_attrs.vsKind='Cartesian'

meshGroup=output_file.create_group('/','meshScaled','')
meshGroup._v_attrs.vsCentering='nodal'
meshGroup._v_attrs.vsIndexOrder='compMinorC'
meshGroup._v_attrs.vsKind='uniform'
#meshGroup._v_attrs.vsLowerBounds=np.array([np.min(FULL_FIELD[1,:,:]),np.min(FULL_FIELD[:,1,:]),np.min(FULL_FIELD[:,:,1])])
meshGroup._v_attrs.vsLowerBounds=mins
meshGroup._v_attrs.vsNumCells=np.array([nX[0]-1,nY[0]-1,nZ2[0]-1])
meshGroup._v_attrs.vsStartCell=np.array([0,0,0])
meshGroup._v_attrs.vsType='mesh'
#meshGroup._v_attrs.vsUpperBounds=np.array([np.max(FULL_FIELD[1,:,:]),np.max(FULL_FIELD[:,1,:]),np.max(FULL_FIELD[:,:,1])])
meshGroup._v_attrs.vsUpperBounds=maxs

timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'

FieldGroup._v_attrs.vsCentering='nodal'
FieldGroup._v_attrs.vsIndexOrder='compMinorC'
FieldGroup._v_attrs.vsLabels='Full_Field'
FieldGroup._v_attrs.vsLimits='globalLimits'
FieldGroup._v_attrs.vsMesh='meshScaled'
FieldGroup._v_attrs.vsTimeGroup='time'
FieldGroup._v_attrs.vsType='variable'

FieldGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
FieldGroup._v_attrs.FXFELSourceFileOrigin='PUFFIN'
FieldGroup._v_attrs.FXFELSourceDirectoryName=os.getcwd()

#Close the file
output_file.close()