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
   print 'Usage: SI_52Puffin <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()


f=tables.open_file(file_name_in,'r')
Electrons=f.root.Particles.read()
n=len(Electrons)
print Electrons.shape

m_X = Electrons[:,0]
m_PX = Electrons[:,1]
m_Y = Electrons[:,2]
m_PY = Electrons[:,3]
m_Z = Electrons[:,4]
m_PZ = Electrons[:,5]
m_WGHT = Electrons[:,6]


#Calculate concentration per/volume
xyz = np.vstack([m_X,m_Y,m_Z]).T

size_x=max(m_X)-min(m_X)
size_y=max(m_Y)-min(m_Y)
size_z=max(m_Z)-min(m_Z)
binnumber=10
cube_volume=(size_x*size_y*size_z)/float(binnumber**3)
print 'Volume of sample cube: ',cube_volume
# print 'Number of particle in record: ',num_of_particles
print 'Size of sample x,y,z: ', size_x,size_y,size_z
print m_WGHT.shape

print (xyz.shape)
#xyz=xyz.T
H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=m_WGHT.flat)
print H.shape

n_p=float(np.amax(H))/cube_volume
print 'Np= ',n_p

# Puffin variables
e_ch=1.602e-19
c=3.0e+8
m=9.11e-31
Pi=3.1415
k_u=251                   # Undulator wave number default=628
a_u=1                     # undulator parameter ? a_u=a_w
e_0=8.854*10**(-12)       # vacuum permitivity
p_tot=np.sqrt((m_PX[:]**2)+(m_PY[:]**2)+(m_PZ[:]**2))

#print np.size(p_tot),p_tot
gamma=(np.sqrt(1+(p_tot/(m*c))**2))
#print 'P_tot/mc ',p_tot/(m*c)
#print np.size(gamma),gamma

omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
gamma_0=np.mean(gamma)
print 'Gamma= ',gamma_0

rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))
lambda_u=(2*Pi)/k_u
lambda_r=(lambda_u/(2*gamma_0**2))*(1+a_u**2)

print 'Rho= ', rho
print 'Lambda_u= ',lambda_u
print 'Lambda_r= ',lambda_r

Lc=lambda_r/(4*Pi*rho)
Lg=lambda_u/(4*Pi*rho)

print 'Lg= ',Lg
print 'Lc= ',Lc
print '4*Pi*Rho= ',4*Pi*rho

# Puffin output arrays:
z2=m_Z/Lc

# Find minimum of z2 and rescale z2 to start from 0.01
min_z2=min(z2)

z2=z2-min_z2+0.01

x_bar=m_X[:]/(np.sqrt(Lg*Lc))
y_bar=m_Y[:]/(np.sqrt(Lg*Lc))
px_bar=m_PX[:]/(m*c*a_u)
py_bar=m_PY[:]/(m*c*a_u)
#sig_gamma_tot
#sig_px_bar
#sig_py_bar
#Ne=m_WGHT[:]/n_p          # weight charge
Ne=m_WGHT[:]          # weight number of electrons

# Combine all read arrays into one
m=np.vstack((x_bar,y_bar,px_bar,py_bar,gamma/gamma_0,z2,Ne)).T

output_file=tables.open_file(file_name_base+'_Puffin.hdf','w')

# Create hdf5 file


# Save the array into hdf5 file
ParticleGroup=output_file.create_array('/','Particles', m)

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
ParticleGroup._v_attrs.vsLabels='x_bar,y_bar,px_bar,py_bar,gamma,z2,Ne'
#Close the file
output_file.close()
out_txt=open(file_name_base+'_Puffin.txt','w')

n=len(m_X)

out_txt.write('__[x_bar]_______[y_bar]_______[px_bar]______[py_bar]_____ [Gamma]________[Z2]__________[Ne]_____'+'\n')
for i in range(n): 
        out_txt.write(      "%.6e" %(x_bar[i]) + \
                            " % .6e" %(y_bar[i]) + \
                            " % .6e" %(px_bar[i]) + \
                            " % .6e" %(py_bar[i]) + \
                            " % .6e" %(gamma[i]/gamma_0) + \
                            " % .6e" %(z2[i]) + \
                            " % .6e" %(Ne[i]) + "\n")
out_txt.close()








