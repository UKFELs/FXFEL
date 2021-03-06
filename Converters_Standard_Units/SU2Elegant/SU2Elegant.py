# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
import numpy as np
import tables
import os
import sys

ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=2.99792458E8

if len(sys.argv)==2:
   file_input_hdf=sys.argv[1]
   print 'Processing file:', file_input_hdf
else:
   print 'Usage: SU2Elegant <input_file> \n'
   sys.exit(1)  
file_name_base  = (file_input_hdf.split('.')[0]).strip()

print 'Loading data from SI HDF5...'

f=tables.open_file(file_input_hdf,'r')
Electrons=f.root.Particles.read()

# Assign particles to arrays and descale them from p/mc to SI

X=Electrons[:,0]
Px=Electrons[:,1]*(me*c)
Y=Electrons[:,2]
Py=Electrons[:,3]*(me*c)
Z=Electrons[:,4]
Pz=Electrons[:,5]*(me*c)
NE=Electrons[:,6]
particleID=np.arange(1,len(X)+1)
print 'Processing data...'

P=np.sqrt(Px**2+Py**2+Pz**2)/(5.34428595e-28*(ERM*1.E6))


Beta_z=Pz/P/(5.34428595e-28*(ERM*1.E6))

xp=Px/Pz
yp=Py/Pz
m_t=Z/(c*Beta_z)
Bunch_Length=np.max(m_t)-np.min(m_t)
Min_Time=np.min(m_t)
m_t_mr=(2.0*Min_Time)+Bunch_Length-m_t

full_array=np.column_stack((X,xp,Y,yp,m_t_mr,P,particleID))

FinalHDF_File=file_name_base+'_Final_HDF.h5'
FinalSDDS_File=file_name_base+'_S2E.sdds'


print 'Creating final HDF5 output...'
output_file=tables.open_file(FinalHDF_File,'w')

ParticleGroup=output_file.create_array('/','spatialPositions',full_array)
del ParticleGroup._v_attrs.CLASS
del ParticleGroup._v_attrs.VERSION
del ParticleGroup._v_attrs.FLAVOR
del ParticleGroup._v_attrs.TITLE
output_file.close()

os.system("hdf2sdds %s %s -dataset=spatialPositions"\
%(FinalHDF_File,FinalSDDS_File))
TotalCharge=sum(NE)*e_ch
print 'Creating SDDS output...'
os.system("sddsconvert %s -rename=columns,D1=x,D2=xp,D3=y,D4=yp,D5=t,D6=p,D7=particleID"\
%(FinalSDDS_File))
os.system("sddsprocess %s %s -define=parameter,Charge,%s,units=C"\
%(FinalSDDS_File,FinalSDDS_File,TotalCharge))
os.system("sddsprocess -convertunits=column,x,m,,1 -convertunits=column,y,m,,1 -convertunits=column,t,s,,1 %s "\
%(FinalSDDS_File))
os.remove(FinalSDDS_File+'~')
os.remove(FinalHDF_File)
print 'Created SDDS: ',FinalSDDS_File