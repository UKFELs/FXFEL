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
import matplotlib.pyplot as plt

import numpy as np
import tables
import sys
import time
import datetime
now = datetime.datetime.now()
print 'Current time: ',now.strftime("%Y-%m-%d %H:%M:%S")

from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.interpolate as interpolate
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

cumulative_nq_Z_FL=np.zeros(0)
cumulative_nq_Y_FL=np.zeros(0)
cumulative_nq_X_FL=np.zeros(0)
cumulative_nq_Z=np.zeros(0)
cumulative_nq_Y=np.zeros(0)
cumulative_nq_X=np.zeros(0)

full_ynew_Z=np.zeros(0)
full_ynew_Y=np.zeros(0)
full_ynew_X=np.zeros(0)


f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()
mA_X = Particles[:,0]
mA_Y = Particles[:,2]    
mA_Z = Particles[:,4]
mA_PX = Particles[:,1]
mA_PY = Particles[:,3]  
mA_PZ = Particles[:,5]
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
DensityFactor=10.0
print 'Total number of electrons: ',TotalNumberOfElectrons
print 'Number of source particles: ',len(mB_X)
print 'Electrons/macroparticles in source data:',round(TotalNumberOfElectrons/len(mB_X))
print 'Desired Electrons/macroparticles in output data:',int(round(TotalNumberOfElectrons/(len(mB_X)*DensityFactor)))

# Set desired particle number per slice - best to set the number that will give 
# an integer when total number of source particles will be divided over it
# Otherwise your particles will be cut at the end
# Very low number can give weird results

NumberOfParticlesPerSlice=100
TotalNumberOfSteps=1+(len(mB_X)/NumberOfParticlesPerSlice)
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
    
    
    #z_axis_length=max(x0_Z)-min(x0_Z)

    #t_knots_z=[min(x0_Z)+0.1*z_axis_length,min(x0_Z)+0.25*z_axis_length,np.mean(x0_Z),max(x0_Z)-0.25*z_axis_length,max(x0_Z)-0.1*z_axis_length]
    #print 'Knots values for LSQ Spline = ',t_knots_z
    #t_knots=[(min(x0_Z)+np.mean(x0_Z))/2,np.mean(x0_Z),(max(x0_Z)+np.mean(x0_Z))/2]
    #print t_knots
    #f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)
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

    full_ynew_Z=np.append(ynew_Z,full_ynew_Z)
    full_ynew_Y=np.append(ynew_Y,full_ynew_Y)
    full_ynew_X=np.append(ynew_X,full_ynew_X)

    #Calculate CDF
cumulative_Z=np.cumsum(full_ynew_Z)
cumulative_Y=np.cumsum(full_ynew_Y)
cumulative_X=np.cumsum(full_ynew_X)
    
    # Remove duplicates from CDF and sort it to avoid chaotic ordering caused by set command
    
    
cumulative_nq_Z_FL=sorted(set(cumulative_Z))
cumulative_nq_Y_FL=sorted(set(cumulative_Y))
cumulative_nq_X_FL=sorted(set(cumulative_X))


    
Num_Of_Target_Particles=int(NumParticlesInSlice*DensityFactor*TotalNumberOfSteps)    
    
print 'Z Cumulative length = ',len(cumulative_nq_Z_FL)
print 'Y Cumulative length = ',len(cumulative_nq_Y_FL)
print 'X Cumulative length = ',len(cumulative_nq_X_FL)
   
xx_0_Z=np.linspace(np.min(mB_Z),np.max(mB_Z),len(cumulative_nq_Z_FL))
xx_0_Y=np.linspace(np.min(mB_Y),np.max(mB_Y),len(cumulative_nq_Y_FL))
xx_0_X=np.linspace(np.min(mB_X),np.max(mB_X),len(cumulative_nq_X_FL))  


cdf_axis_length=max(cumulative_nq_Z_FL)-min(cumulative_nq_Z_FL)

tt_knots_z=[min(cumulative_nq_Z_FL)+0.1*cdf_axis_length,min(cumulative_nq_Z_FL)+0.25*cdf_axis_length,np.mean(cumulative_nq_Z_FL),max(cumulative_nq_Z_FL)-0.25*cdf_axis_length,max(cumulative_nq_Z_FL)-0.1*cdf_axis_length]
#tt_knots_z=[min(cumulative_nq_Z_FL)+0.1*cdf_axis_length,np.mean(cumulative_nq_Z_FL),max(cumulative_nq_Z_FL)-0.1*cdf_axis_length]
print 'Knots values for LSQ Spline in CDF = ',tt_knots_z

ff_Z = interpolate.LSQUnivariateSpline(cumulative_nq_Z_FL, xx_0_Z,tt_knots_z)
    
#ff_Z = interpolate.UnivariateSpline(cumulative_nq_Z_FL, xx_0_Z)
ff_Y = InterpolatedUnivariateSpline(cumulative_nq_Y_FL, xx_0_Y)
ff_X = InterpolatedUnivariateSpline(cumulative_nq_X_FL, xx_0_X)

plt.title('Cumulative density function')
plt.xlabel('Z')
plt.ylabel('CDF')
plt.grid(True)
plt.plot(xx_0_Z, cumulative_nq_Z_FL, label='CDF')
plt.show()

 
#density_Z=np.linspace(min(cumulative_nq_Z_FL),max(cumulative_nq_Z_FL),Num_Of_Target_Particles)
density_Z=np.random.uniform(low=min(cumulative_nq_Z_FL), high=max(cumulative_nq_Z_FL), size=(Num_Of_Target_Particles))
density_Y=np.random.uniform(low=min(cumulative_nq_Y_FL), high=max(cumulative_nq_Y_FL), size=(Num_Of_Target_Particles))
density_X=np.random.uniform(low=min(cumulative_nq_X_FL), high=max(cumulative_nq_X_FL), size=(Num_Of_Target_Particles))

print 'Target number of particles = ',Num_Of_Target_Particles

Full_Z=ff_Z(density_Z)
Full_Y=ff_Y(density_Y)
Full_X=ff_X(density_X)
print 'New Z length = ',len(Full_Z)    
print 'New Y length = ',len(Full_Y)
print 'New X length = ',len(Full_X)

    # Merge particles to common array (all steps together)


print 'Interpolating momentum data on new grid. This may take a while...'

Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')

No_Particles_Per_Record=np.zeros((len(Full_Z),1))

No_Particles_Per_Record[0:len(Full_Z)]=int(round((TotalNumberOfElectrons/len(Full_Z))))


x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,No_Particles_Per_Record.flat]).T

#x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_Y.flat,Full_Z.flat,No_Particles_Per_Record.flat]).T



plt.title('Density profiles')
plt.xlabel('Z')
plt.ylabel('Density')
plt.hist(Full_Z,weights=No_Particles_Per_Record,bins=100)
#plt.plot(x0_Z,y0_Z,'yo')
plt.ticklabel_format(style='sci')
plt.grid(True)
#plt.plot(xnew_Z_plt,ynew_Z_plt,linewidth=4, color='red')
plt.show()

output_file=tables.open_file(file_name_base+'_TST.si5','w')

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
#ParticleGroup._v_attrs.vsLabels='x,y,z,NE'
ParticleGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
ParticleGroup._v_attrs.FXFELDensityFactor=DensityFactor
#Close the file
output_file.close()

