# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:03:14 2015

@author: piotrt
"""
# Import necessary libraries
import tables

import numpy as np
import gc
import sys


if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SU2Puffin <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


f=tables.open_file(file_name_in,'r')
Electrons=f.root.Particles.read()
n=len(Electrons)

c=3.0e+8
m=9.11e-31
e_ch=1.602e-19              # charge of one electron

# Assign the data to array and descale the units from p/mc to SI !

m_X = Electrons[:,0]
m_PX = Electrons[:,1]*(m*c)
m_Y = Electrons[:,2]
m_PY = Electrons[:,3]*(m*c)
m_Z = Electrons[:,4]
m_PZ = Electrons[:,5]*(m*c)
m_WGHT = Electrons[:,6]

del Electrons
gc.collect()


size_x=max(m_X)-min(m_X)
size_y=max(m_Y)-min(m_Y)
size_z=max(m_Z)-min(m_Z)
binnumber=20
# print 'Number of particle in record: ',num_of_particles
#print 'Size of sample x,y,z: ', size_x,size_y,size_z
#print m_WGHT.shape

Hn, edgesn = np.histogram(m_Z, bins = binnumber,normed=False,weights=m_WGHT.flat)
#print Hn
if all(x==Hn[0] for x in Hn)==True:
    print 'This seems like a Flat-Top beam...'
    z_low  = np.max(m_Z)-(size_z*float((binnumber/2))+1.0)
    z_high = np.max(m_Z)-(size_z*float((binnumber/2))-1.0)
else:
    print 'Not a Flat-Top beam...'
    z_low  =np.max(m_Z)-(size_z*float(np.argwhere(Hn == Hn.max()))/(binnumber-1.0))
    z_high =np.max(m_Z)-(size_z*float(np.argwhere(Hn == Hn.max()))/(binnumber+1.0))
#print 'Max value for x = ', np.argwhere(Hn == Hn.max())
#print 'Edges = ',edgesn[np.argwhere(Hn == Hn.max())]
#print edgesn

# Select slice from the beam with the peak number of electrons


#z_low  =np.max(m_Z)-(size_z*float(np.argwhere(Hn == Hn.max()))/(binnumber-1.0))
#z_high =np.max(m_Z)-(size_z*float(np.argwhere(Hn == Hn.max()))/(binnumber+1.0))
#print 'Zlow = ',z_low,' Zhigh = ',z_high
mA_X=m_X[(m_Z>=z_low) & (m_Z<=z_high)]
mA_Y=m_Y[(m_Z>=z_low) & (m_Z<=z_high)]

from scipy.spatial import ConvexHull
points = np.vstack([mA_X,mA_Y]).T
#print points   # 30 random points in 2-D
hull = ConvexHull(points)
#print 'Area = ',hull.volume


#==============================================================================
# # Plot the selected convex hull
# import matplotlib.pyplot as plt
# plt.plot(points[:,0], points[:,1], 'o')
# for simplex in hull.simplices:
#      plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
# plt.show()
#==============================================================================

# Calculate peak current
Peak_Current=((np.max(Hn)*e_ch)/(size_z/binnumber))*c
print 'Peak current = ',Peak_Current
Peak_Concentration=Peak_Current/c/e_ch/(hull.volume)
print 'Peak concentration = ',Peak_Concentration

gc.collect()

# Puffin variables
e_ch=1.602e-19
c=299792458.0
m=9.11e-31
Pi=np.pi
k_u=228.4727               # Undulator wave number default=628 k_u=2*Pi/l_w
a_u=1.0121809              # undulator parameter ? a_u=a_w
e_0=8.854E-12              # vacuum permitivity
p_tot=np.sqrt((m_PX[:]**2.0)+(m_PY[:]**2.0)+(m_PZ[:]**2.0))

#print np.size(p_tot),p_tot
gamma=(np.sqrt(1.0+(p_tot/(m*c))**2))
#print 'P_tot/mc ',p_tot/(m*c)
#print np.size(gamma),gamma

omega_p=np.sqrt((e_ch*e_ch*Peak_Concentration)/(e_0*m))
gamma_0=np.mean(gamma)
print 'Gamma= ',gamma_0


rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))
#Temporary change of rho to rho = 0.005
#rho = 0.0050

lambda_u=(2*Pi)/k_u

# Calculated for planar undulator -> (1+(a^2)/2)
lambda_r=(lambda_u/(2.0*gamma_0**2.0))*(1+(a_u**2.0)/2.0)
#lambda_r=(lambda_u/(2*gamma_0**2))*(1+a_u**2)

print 'Rho= ', rho
print 'Lambda_u= ',lambda_u
print 'Lambda_r= ',lambda_r

Lc=lambda_r/(4*Pi*rho)
Lg=lambda_u/(4*Pi*rho)

print 'Lg= ',Lg
print 'Lc= ',Lc
#print '4*Pi*Rho= ',4*Pi*rho

z2=m_Z/Lc

# Find minimum of z2 and rescale z2 to start from 0.01
# Flip the beam to keep proper direction in Puffin (z-ct)
min_z2=min(z2)
mean_z2=np.mean(z2)
z2=(mean_z2-z2)+min_z2
min_z22=min(z2)
z2=z2-min_z22+0.01

x_bar=m_X[:]/(np.sqrt(Lg*Lc))
y_bar=m_Y[:]/(np.sqrt(Lg*Lc))
px_bar=m_PX[:]/(m*c*a_u)
py_bar=-1.0*m_PY[:]/(m*c*a_u)


# Centering the beam - uncomment if necessary
#avg_x=np.mean(x_barN)
#avg_y=np.mean(y_barN)

#x_bar=avg_x-x_barN
#y_bar=avg_y-y_barN

Ne=m_WGHT[:]/(Peak_Concentration*Lg*Lc*Lc)          # weight number of electrons and scale with Lg*Lc^2

# Scale electrons weights wit scaled_n_peak (given rho in use !!)
#Ne=Ne/scaled_n_peak
del m_X, m_PX, m_Y, m_PY, m_Z, m_PZ, m_WGHT,p_tot 
gc.collect()

# Combine all read arrays into one
m_Arr=np.vstack((x_bar,y_bar,z2,px_bar,py_bar,gamma/gamma_0,Ne)).T

del x_bar,y_bar,px_bar,py_bar,gamma,z2,Ne
gc.collect()

output_file=tables.open_file(file_name_base+'_Puffin.hdf','w')
# Create hdf5 file

print 'Max Z2 value = ',np.max(m_Arr[:,2])
# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','electrons', m_Arr)

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
ParticleGroup._v_attrs.vsLabels='x_bar,y_bar,z2,px_bar,py_bar,gamma,Ne'
#Close the file
output_file.close()








