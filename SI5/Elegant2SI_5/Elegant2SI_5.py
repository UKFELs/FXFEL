# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
import numpy as np
import os
import csv
import sys
import tables
import datetime
now = datetime.datetime.now()
print 'Conversion time: ',now.strftime("%Y-%m-%d %H:%M:%S")

ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=3e+08

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Elegant2SI_5 <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()




TotalCharge=os.popen("sdds2stream -parameters=Charge %s" %(file_name_in)).read()

print 'Total Charge= ',float(TotalCharge)
os.system('sdds2stream -column=x,xp,y,yp,t,p %s > %s'%(file_name_in,file_name_base+'OUT'))

linecount = sum(1 for line in open(file_name_base+'OUT'))

ParticleChargePerRecord=float(TotalCharge)/linecount
print 'Single Particle Charge= ',float(ParticleChargePerRecord)

m_X = np.zeros((linecount,1), dtype=np.float64)
m_XP = np.zeros((linecount,1), dtype=np.float64)
m_Y = np.zeros((linecount,1), dtype=np.float64)
m_YP = np.zeros((linecount,1), dtype=np.float64)
m_T = np.zeros((linecount,1), dtype=np.float64)
m_P = np.zeros((linecount,1), dtype=np.float64)


f = csv.reader(open(file_name_base+'OUT', 'r'),delimiter=' ')
for i, row in enumerate(f):
 #   m[i] = float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6])
     m_X[i] = float(row[0]),
     m_XP[i] = float(row[1]),
     m_Y[i] = float(row[2]),
     m_YP[i] = float(row[3]),
     m_T[i] = float(row[4]),
     m_P[i] = float(row[5]),




m_PZ=(5.36E-28*(ERM*1.E6)*m_P)/(np.sqrt((m_XP**2)+(m_YP**2)+1))
m_PX=m_PZ*m_XP
m_PY=m_PZ*m_YP

#m_ZP=np.sqrt((m_P[:]**2)-(m_XP[:]**2)-(m_YP[:]**2))
#Beta_z=np.sqrt((m_ZP**2+m_XP**2+m_YP**2)/m_ZP**2)

Beta_z=m_PZ/(np.sqrt(m_PX**2+m_PY**2+m_PZ**2))

p_tot=np.sqrt((m_PX[:]**2)+(m_PY[:]**2)+(m_PZ[:]**2))
Gamma=(np.sqrt(1+(p_tot/(me*c))**2))
V_z=m_PZ/(Gamma*me)
Z_true=V_z*m_T


# By some reason data is flipped along Z-axis, to correct it we need to flip it again,
# therefore below lines that find min Z values and rotate the beam in this point (mirror flip).
# So e.g. beam instead of looking like <<o will look like o>>

m_Z=m_T*c*Beta_z
mirror_flip_z=np.min(m_Z)
m_Z_mr=m_Z-((m_Z-mirror_flip_z)*2)

No_Particles_Per_Record=np.zeros((linecount,1))

No_Particles_Per_Record[0:len(m_Z)]=(ParticleChargePerRecord/e_ch)

x_px_y_py_z_pz_NE = np.hstack([m_X,m_PX,m_Y,m_PY,m_Z_mr,m_PZ,No_Particles_Per_Record])




output_file=tables.open_file(file_name_base+'.si5','w')

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
ParticleGroup._v_attrs.FXFELSourceFileOrigin='ELEGANT'
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in

#Close the file
output_file.close()

os.remove(file_name_base+'OUT')


