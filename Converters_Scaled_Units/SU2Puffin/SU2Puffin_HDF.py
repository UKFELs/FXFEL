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
import SUF  # SU Format


def SU2Puffin(fnamein):

    file_name_base  = (fnamein.split('.')[0]).strip()

    m_X, m_PX, m_Y, m_PY, m_Z, m_PZ, m_WGHT = SUF.readSUF(fnamein)

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

    del xyz
    gc.collect()

# Puffin variables
    e_ch=1.602e-19
    c=3.0e+8
    m=9.11e-31
    c0 = 2.99792458e8
    qe = 1.60217653e-19
    eps0 = 8.854187817e-12
    me = 9.1093826e-31
    h = 6.626e-34
    Pi=3.1415
    k_u=228.47946               # Undulator wave number default=628 k_u=2*Pi/l_w
    a_u=0.71572                     # undulator parameter ? a_u=a_w
    e_0=8.854E-12              # vacuum permitivity
    
    
    p_tot=np.sqrt((m_PX[:]**2)+(m_PY[:]**2)+(m_PZ[:]**2))

    gamma=(np.sqrt(1+(p_tot)**2))

    omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
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
    print '4*Pi*Rho= ',4*Pi*rho

# The below line allow user to force different value of rho than calculated
# for purpose of electrons_weight scaling in HDF files in Puffin.
# given_rho=0.005
#given_rho=rho
#n_peak=(e_0*m/(e_ch**2.0))*(((given_rho*gamma_0)**(3.0/2.0)*(4.0*c*k_u))/a_u)**2.0
    print 'Peak density = ',n_p
#print 'Calculated n_peak from rho = ',n_peak
#scaled_n_peak=n_peak*Lg*Lc**2.0
#print 'Scaled n_peak = ',scaled_n_peak
# Puffin output arrays:
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
    px_bar=m_PX[:]/a_u
    py_bar=-1.0*m_PY[:]/a_u


# Centering the beam - uncomment if necessary
#avg_x=np.mean(x_barN)
#avg_y=np.mean(y_barN)

#x_bar=avg_x-x_barN
#y_bar=avg_y-y_barN



#sig_gamma_tot
#sig_px_bar
#sig_py_bar
#Ne=m_WGHT[:]/n_p          # weight charge
    Ne=m_WGHT[:]/(n_p*Lg*Lc*Lc)          # weight number of electrons and scale with Lg*Lc^2

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

# Close the file
    output_file.close()


if __name__ == '__main__':

    if len(sys.argv)==2:
        fname = sys.argv[1]
        print 'Processing file:', fname
        SU2Puffin(fname)
    else:
        print 'Usage: SU2Puffin <FileName> \n'
        sys.exit(1)
        
        

