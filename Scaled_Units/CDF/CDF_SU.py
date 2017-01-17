# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:07:49 2016

@author: piotrt
"""
# This scripts uses Cumulative Distribution Function to increase number of particles
# used in FEL simulations.
# The first what has to be done is to sort particles along Z-axis, the purpose of doing so is 
# to easy access the ordered particles in loop 'walking' along the Z-axis. This method saves us
# from reaching low particle number in slice what could happen if we just do the slicing
# along Z-axis using length values (splitting the beam to smaller parts of certain length). Each slice
# has same number of particles.
# The next step is to generate density histogram - this is done by using 'histogramdd' for
# each axis separatyly and also using particle number per macroparticle as weights - I assume
# that sometimes source file will have macroparticle of where electron numbers per particle is not the same
# Then the shape is approximated by linear function and next, proper CDF function is created
# Finally, the CDF function is projected on randomly placed x,y,z values - the Z-axis values are supposed to
# be uniformely spaced that is why 'linspace' is used and not 'random.uniform'.
# The above is done in loop for each slice. The results are appended to larger common array.
# The last step for generating the particle beam is to project momentum onto new set of particles.
# This is done using 'interpolate.griddata' what projects exact shape of momentum. The number of
# electrons per macroparticles is equalized - so the total sum of electrons is same as at the beginnig
# just divided over larger number of macroparticles.
# The output file is saved as HDF5 file with VizSchema metadata added.


import numpy as np
import tables
import sys
import time
import datetime
now = datetime.datetime.now()
print 'Current time: ',now.strftime("%Y-%m-%d %H:%M:%S")

from scipy.interpolate import InterpolatedUnivariateSpline

e_ch=1.602e-19


if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: DensityCalc <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()

 
# time - just for benchmarking, not necessary
start = time.time()
def elapsed():
    return time.time() - start


# pre-allocate array and load data into array
Full_X=np.zeros(0)
Full_PX=np.zeros(0)
Full_Y=np.zeros(0)
Full_PY=np.zeros(0)
Full_Z=np.zeros(0)
Full_PZ=np.zeros(0)

c=3.0e+8                    # Speed of light
m=9.11e-31                  # mass of electron



f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()
mA_X = Particles[:,0]
mA_Y = Particles[:,2]    
mA_Z = Particles[:,4]
# Descale the momentum units from p/mc to SI - to keep compatibility with rest of calculations
mA_PX = Particles[:,1]*(m*c)
mA_PY = Particles[:,3]*(m*c)  
mA_PZ = Particles[:,5]*(m*c)
mA_WGHT = Particles[:,6]


xyzW = np.vstack((mA_X.flat,mA_Y.flat,mA_Z.flat,mA_WGHT.flat)).T

print np.shape(xyzW)
xyzW=xyzW[xyzW[:,2].argsort()]

# Just the time for benchmark
print '%.3fs: loaded' % elapsed()
#Calculate total charge
TotalNumberOfElectrons=np.sum(mA_WGHT)



mB_X=xyzW[:,0].flat
mB_Y=xyzW[:,1].flat
mB_Z=xyzW[:,2].flat
mB_WGHT=xyzW[:,3].flat
print len(mB_X)


# Set the multiplier of desired particle number
# The highher the number the more time it will take
# Please not that if you set this too high the algorithm will still work
# but at the end your particles will have smaller than one what means that 
# you have less than one electron per particle
# Currently there is no safety mechanism to avoid this
DensityFactor=1.0
print 'Total number of electrons: ',TotalNumberOfElectrons
print 'Number of source particles: ',len(mB_X)
print 'Electrons/macroparticles in source data:',round(TotalNumberOfElectrons/len(mB_X))
print 'Desired Electrons/macroparticles in output data:',int(round(TotalNumberOfElectrons/(len(mB_X)*DensityFactor)))

# Set desired particle number per slice - best to set the number that will give 
# an integer when total number of source particles will be divided over it
# Otherwise your particles will be cut at the end
# Very low number can give weird results

NumberOfParticlesPerSlice=50
TotalNumberOfSteps=len(mB_X)/NumberOfParticlesPerSlice
#part_no_in_step=len(mB_X)/total_no_steps

print 'Total number of steps= ',TotalNumberOfSteps
print 'Macroparticles in step= ',NumberOfParticlesPerSlice



for stepnumber in range(1,TotalNumberOfSteps):
    m_X=mB_X[(stepnumber-1)*NumberOfParticlesPerSlice:NumberOfParticlesPerSlice*stepnumber]
#    print len(m_X)
    m_Y=mB_Y[(stepnumber-1)*NumberOfParticlesPerSlice:NumberOfParticlesPerSlice*stepnumber]
    m_Z=mB_Z[(stepnumber-1)*NumberOfParticlesPerSlice:NumberOfParticlesPerSlice*stepnumber]
    m_WGHT=mB_WGHT[(stepnumber-1)*NumberOfParticlesPerSlice:NumberOfParticlesPerSlice*stepnumber]
#    print 'Stepnumber= ',stepnumber      
    NumParticlesInSlice=len(m_X)


   
    binnumber=10


    Hz, edges_Z = np.histogramdd(m_Z, bins = binnumber,normed=False,weights=m_WGHT.flat)
    Hy, edges_Y = np.histogramdd(m_Y, bins = binnumber,normed=False,weights=m_WGHT.flat)
    Hx, edges_X = np.histogramdd(m_X, bins = binnumber,normed=False,weights=m_WGHT.flat)


    y0_Z = Hz
    y0_Y = Hy
    y0_X = Hx


    # Linear approximation of density function
    x0_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber)
    x0_Y = np.linspace(np.min(m_Y),np.max(m_Y),binnumber)
    x0_X = np.linspace(np.min(m_X),np.max(m_X),binnumber)

    from scipy.interpolate import interp1d
    f_Z = interp1d(x0_Z, y0_Z)
    f_Y = interp1d(x0_Y, y0_Y)
    f_X = interp1d(x0_X, y0_X)


    xnew_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber)
    xnew_Y = np.linspace(np.min(m_Y),np.max(m_Y),binnumber)
    xnew_X = np.linspace(np.min(m_X),np.max(m_X),binnumber)



    ynew_Z=f_Z(xnew_Z)
    ynew_Y=f_Y(xnew_Y)
    ynew_X=f_X(xnew_X)


    ynew_Z=ynew_Z.clip(min=0)
    ynew_Y=ynew_Y.clip(min=0)
    ynew_X=ynew_X.clip(min=0)

    #Scale values of function so, that the area under function is 1

    ynew_Z=ynew_Z/np.sum(ynew_Z)
    ynew_Y=ynew_Y/np.sum(ynew_Y)
    ynew_X=ynew_X/np.sum(ynew_X)

    #Calculate CDF

    cumulative_Z=np.cumsum(ynew_Z)
    cumulative_Y=np.cumsum(ynew_Y)
    cumulative_X=np.cumsum(ynew_X)
    
    # Remove duplicates from CDF and sort it to avoid chaotic ordering caused by set command

    cumulative_nq_Z=sorted(set(cumulative_Z))
    cumulative_nq_Y=sorted(set(cumulative_Y))
    cumulative_nq_X=sorted(set(cumulative_X))



    xx_0_Z=np.linspace(np.min(m_Z),np.max(m_Z),len(cumulative_nq_Z))
    xx_0_Y=np.linspace(np.min(m_Y),np.max(m_Y),len(cumulative_nq_Y))
    xx_0_X=np.linspace(np.min(m_X),np.max(m_X),len(cumulative_nq_X))

    ff_Z = interp1d(cumulative_nq_Z, xx_0_Z)
    ff_Y = interp1d(cumulative_nq_Y, xx_0_Y)
    ff_X = interp1d(cumulative_nq_X, xx_0_X)

    Num_Of_Target_Particles=int(NumParticlesInSlice*DensityFactor)
 
    density_Z=np.linspace(min(cumulative_nq_Z),max(cumulative_nq_Z),Num_Of_Target_Particles)
    density_Y=np.random.uniform(low=min(cumulative_nq_Y), high=max(cumulative_nq_Y), size=(Num_Of_Target_Particles))
    density_X=np.random.uniform(low=min(cumulative_nq_X), high=max(cumulative_nq_X), size=(Num_Of_Target_Particles))

    z_density_Z=ff_Z(density_Z)
    z_density_Y=ff_Y(density_Y)
    z_density_X=ff_X(density_X)
    


    # Merge particles to common array (all steps together)
    Full_X=np.append(z_density_X,Full_X)   
    Full_Y=np.append(z_density_Y,Full_Y)
    Full_Z=np.append(z_density_Z,Full_Z)

print 'Interpolating momentum data on new grid. This may take a while...'
import scipy.interpolate as interpolate
Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')

No_Particles_Per_Record=np.zeros((len(Full_Z),1))

No_Particles_Per_Record[0:len(Full_Z)]=int(round((TotalNumberOfElectrons/len(Full_Z))))


# Scale units from SI to p = p/mc

Full_PX=Full_PX/(m*c)
Full_PY=Full_PY/(m*c)
Full_PZ=Full_PZ/(m*c)

x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,No_Particles_Per_Record.flat]).T




output_file=tables.open_file(file_name_base+'_Dense.si5','w')

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
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
ParticleGroup._v_attrs.FXFELDensityFactor=DensityFactor
#Close the file
output_file.close()

