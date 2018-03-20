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
import getTwiss
import matchTwiss
from puffDataClass import puffData
from undulator import undulator

def SU2Matched(fnamein):

#    Store basename of file

    file_name_base  = (fnamein.split('.')[0]).strip()

    puffVars = puffData()

# Initialize CLARA base parameters

    puffVars.aw = 0.8745*np.sqrt(2.)   # The PEAK!!!
    puffVars.gamma0 = 489.237
    puffVars.lw = 0.025
    puffVars.rho = 0.005
    puffVars.undtype = 'planepole'
    puffVars.ux = 0.
    puffVars.uy = 1.

    emitx = 1.022e-9
    emity = 1.022e-9

# Generate the rest of the Puffin scaling from the above

    puffVars.genParams()  # generate rest of scaled params

#    Get beam distribution and scale

    x, px, y, py, z, pz, wghts = SUF.readSUF(fnamein)

    p_tot=np.sqrt((px[:]**2)+(py[:]**2)+(pz[:]**2))


    gamma = (np.sqrt(1+(p_tot)**2))
    gamma = gamma / puffVars.gamma0

    xb = x / (np.sqrt(puffVars.lg * puffVars.lc))
    yb = y / (np.sqrt(puffVars.lg * puffVars.lc))

    pxb = px[:]/(puffVars.me * puffVars.c0 * puffVars.au)
    pyb = -1.0 * py[:]/(puffVars.me * puffVars.c0 * puffVars.au)

    z2 = -z / puffVars.lc   # (ct - z) / lc
    z2 = z2 - min(z2) + 0.001

    wghts = wghts / puffVars.npkbar


    oname = file_name_base+'_matched.h5'

    MPs=np.vstack((xb, pxb, yb, pyb, z2, gamma, wghts)).T

    output_file=tables.open_file(oname,'w')


    ParticleGroup=output_file.create_array('/','Particles', MPs)
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
        
        

