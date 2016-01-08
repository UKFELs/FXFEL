# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
import numpy as np
import os
import csv
import sys

ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=3e+08

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Ele2Masp <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()




TotalCharge=os.popen("sdds2stream -parameters=Charge %s" %(file_name_in)).read()

print 'Total Charge= ',float(TotalCharge)
os.system('sdds2stream -column=x,xp,y,yp,t,p %s > %s'%(file_name_in,file_name_base+'OUT'))

linecount = sum(1 for line in open(file_name_base+'OUT'))

ParticleChargePerRecord=float(TotalCharge)/linecount
print 'Single Particle Charge= ',float(ParticleChargePerRecord)

m_X = np.zeros((linecount,1), dtype=np.float32)
m_XP = np.zeros((linecount,1), dtype=np.float32)
m_Y = np.zeros((linecount,1), dtype=np.float32)
m_YP = np.zeros((linecount,1), dtype=np.float32)
m_T = np.zeros((linecount,1), dtype=np.float32)
m_P = np.zeros((linecount,1), dtype=np.float32)


f = csv.reader(open(file_name_base+'OUT', 'r'),delimiter=' ')
for i, row in enumerate(f):
 #   m[i] = float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6])
     m_X[i] = float(row[0]),
     m_XP[i] = float(row[1]),
     m_Y[i] = float(row[2]),
     m_YP[i] = float(row[3]),
     m_T[i] = float(row[4]),
     m_P[i] = float(row[5]),




m_PZ=(m_P/(((np.sqrt((m_XP**2)+(m_YP**2)+1))))*(ERM*1.E6))*5.36*10**(-28)
m_PX=m_PZ*m_XP*5.36*10**(-28)
m_PY=m_PZ*m_YP*5.36*10**(-28)

#m_ZP=np.sqrt((m_P[:]**2)-(m_XP[:]**2)-(m_YP[:]**2))
#Beta_z=np.sqrt((m_ZP**2+m_XP**2+m_YP**2)/m_ZP**2)

Beta_z=m_PZ/(np.sqrt(m_PX**2+m_PY**2+m_PZ**2))
m_Z=m_T*c*Beta_z

out_txt=open(file_name_base+'_ERM.txt','w')

for i in range(linecount): 
        out_txt.write(      "% .4E" %(m_X[i]) + \
                            " % .4E" %(m_PX[i]) + \
                            " % .4E" %(m_Y[i]) + \
                            " % .4E" %(m_PY[i]) + \
                            " % .4E" %(m_Z[i]) + \
                            " % .4E" %(m_PZ[i]) + \
                            " % .4E" %((ParticleChargePerRecord)/e_ch) + "\n")
out_txt.close()
os.remove(file_name_base+'OUT')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(m_X, m_Y, m_Z,c=m_PZ)
plt.show()

