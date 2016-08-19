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
import matplotlib.pyplot as plt
now = datetime.datetime.now()
print 'Current time: ',now.strftime("%Y-%m-%d %H:%M:%S")


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



f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()
mA_X = Particles[:,0]
mA_Y = Particles[:,2]    
mA_Z = Particles[:,4]
mA_PX = Particles[:,1]
mA_PY = Particles[:,3]  
mA_PZ = Particles[:,5]
mA_WGHT = Particles[:,6]



#*************************************************************
# The below section calculates size of the bin according
# to value of Lc (binnumbers=total_length/Lc)

Pi=3.1415
k_u=157.075
a_u=1
c=3.0e+8
m=9.11e-31
e_0=8.854E-12 

xyz = np.vstack([mA_X,mA_Y,mA_Z]).T
size_x=max(mA_X)-min(mA_X)
size_y=max(mA_Y)-min(mA_Y)
size_z=max(mA_Z)-min(mA_Z)
binnumber=100
cube_volume=(size_x*size_y*size_z)/float(binnumber**3)
print 'Volume of sample cube: ',cube_volume
# print 'Number of particle in record: ',num_of_particles
print 'Size of sample x,y,z: ', size_x,size_y,size_z

H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=mA_WGHT.flat)
print H.shape

n_p=float(np.amax(H))/cube_volume
print 'Np= ',n_p


p_tot=np.sqrt((mA_PX[:]**2)+(mA_PY[:]**2)+(mA_PZ[:]**2))
gamma=(np.sqrt(1+(p_tot/(m*c))**2))
gamma_0=np.mean(gamma)
omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))
lambda_u=(2*Pi)/k_u
lambda_r=(lambda_u/(2*gamma_0**2))*(1+a_u**2)
Lc=lambda_r/(4*Pi*rho)
number_of_bins=int((size_z/Lc)*1.10)
print 'Lc = ',Lc
print 'Number of bins = ',number_of_bins

# End of bin size calculations
#*************************************************************

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


# Set the multiplier of desired particle number
# The highher the number the more time it will take
# Please not that if you set this too high the algorithm will still work
# but at the end your particles will have smaller than one what means that 
# you have less than one electron per particle
# Currently there is no safety mechanism to avoid this
DensityFactor=10.0
print 'Initial charge of particles = ',TotalNumberOfElectrons*e_ch
print 'Total number of electrons: ',TotalNumberOfElectrons
print 'Number of source macroparticles: ',len(mB_X)
print 'Electrons/macroparticles in source data:',round(TotalNumberOfElectrons/len(mB_X))
print 'Desired Electrons/macroparticles in output data:',int(round(TotalNumberOfElectrons/(len(mB_X)*DensityFactor)))

# Set desired particle number per slice - best to set the number that will give 
# an integer when total number of source particles will be divided over it
# Otherwise your particles will be cut at the end
# Very low number can give weird results

NumberOfParticlesPerSlice=len(mB_X)


print 'Macroparticles in job= ',NumberOfParticlesPerSlice



 
m_X=mB_X[:]
m_Y=mB_Y[:]
m_Z=mB_Z[:]
m_WGHT=mB_WGHT[:]
          


binnumber_Z=number_of_bins   
#binnumber_Z=100
binnumber_X=100
binnumber_Y=100

print'Binnumber X,Y,Z = ',binnumber_X,binnumber_Y,binnumber_Z
Hz, edges_Z = np.histogramdd(m_Z, bins = binnumber_Z,normed=False,weights=m_WGHT.flat)
Hy, edges_Y = np.histogramdd(m_Y, bins = binnumber_Y,normed=False,weights=m_WGHT.flat)
Hx, edges_X = np.histogramdd(m_X, bins = binnumber_X,normed=False,weights=m_WGHT.flat)


y0_Z = Hz
y0_Y = Hy
y0_X = Hx


    # Linear approximation of density function
x0_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber_Z)
x0_Y = np.linspace(np.min(m_Y),np.max(m_Y),binnumber_Y)
x0_X = np.linspace(np.min(m_X),np.max(m_X),binnumber_X)
    
from scipy import interpolate
from scipy.interpolate import interp1d
x_axis_length=max(x0_X)-min(x0_X)
y_axis_length=max(x0_Y)-min(x0_Y)
z_axis_length=max(x0_Z)-min(x0_Z)

t_knots_x=[min(x0_X)+0.1*x_axis_length,min(x0_X)+0.25*x_axis_length,np.mean(x0_X),max(x0_X)-0.25*x_axis_length,max(x0_X)-0.1*x_axis_length]
t_knots_y=[min(x0_Y)+0.1*y_axis_length,min(x0_Y)+0.25*y_axis_length,np.mean(x0_Y),max(x0_Y)-0.25*y_axis_length,max(x0_Y)-0.1*y_axis_length]
t_knots_z=[min(x0_Z)+0.1*z_axis_length,min(x0_Z)+0.25*z_axis_length,np.mean(x0_Z),max(x0_Z)-0.25*z_axis_length,max(x0_Z)-0.1*z_axis_length]
#t_knots_z=[min(x0_Z)+0.25*z_axis_length,np.mean(x0_Z),max(x0_Z)-0.25*z_axis_length]
#t_knots_z=[(min(x0_Z)+np.mean(x0_Z))/2,np.mean(x0_Z),(max(x0_Z)+np.mean(x0_Z))/2]
print 'Knots values for LSQ Spline = ',t_knots_z

f_X = interpolate.LSQUnivariateSpline(x0_X, y0_X,t_knots_x)
f_Y = interpolate.LSQUnivariateSpline(x0_Y, y0_Y,t_knots_y)
f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)
#f_Z = interpolate.UnivariateSpline(x0_Z, y0_Z)
#f_Z = interp1d(x0_Z, y0_Z)
#f_Y = interp1d(x0_Y, y0_Y)
#f_X = interp1d(x0_X, y0_X)





xnew_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber_Z)
#xnew_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber_Z*10)
xnew_Y = np.linspace(np.min(m_Y),np.max(m_Y),binnumber_Y)
xnew_X = np.linspace(np.min(m_X),np.max(m_X),binnumber_X)



ynew_Z=f_Z(xnew_Z)
ynew_Y=f_Y(xnew_Y)
ynew_X=f_X(xnew_X)

xnew_Z_plt=xnew_Z[:]
ynew_Z_plt=ynew_Z[:]


my_dpi=96
plt.figure(figsize=(1500/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.title('Density profile')
plt.xlabel('Z')
plt.ylabel('Density')
plt.plot(x0_Z,y0_Z,'yo',label='Initial data density values')
plt.ticklabel_format(style='sci')
plt.grid(True)
plt.plot(xnew_Z,ynew_Z,linewidth=4, color='red',label='Fitted density profile')
plt.legend(loc='upper right')
plt.savefig('Initial_density_plot.png',dpi=300)
plt.show()

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


my_dpi=96
plt.figure(figsize=(1500/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.title('Cumulative density function')
plt.xlabel('Z')
plt.ylabel('CDF')
plt.grid(True)
plt.plot(xx_0_Z, cumulative_nq_Z, label='Calculated CDF')
plt.legend(loc='upper right')
plt.savefig('CDF_plot.png',dpi=300)
plt.show()




#cdf_axis_length_X=max(cumulative_nq_X)-min(cumulative_nq_X)
#cdf_axis_length_Y=max(cumulative_nq_Y)-min(cumulative_nq_Y)
cdf_axis_length_Z=max(cumulative_nq_Z)-min(cumulative_nq_Z)

#tt_knots_x=[min(cumulative_nq_X)+0.1*cdf_axis_length_X,min(cumulative_nq_X)+0.25*cdf_axis_length_X,np.mean(cumulative_nq_X),\
#max(cumulative_nq_X)-0.25*cdf_axis_length_X,max(cumulative_nq_X)-0.1*cdf_axis_length_X]

#tt_knots_y=[min(cumulative_nq_Y)+0.1*cdf_axis_length_Y,min(cumulative_nq_Y)+0.25*cdf_axis_length_Y,np.mean(cumulative_nq_Y),\
#max(cumulative_nq_Y)-0.25*cdf_axis_length_Y,max(cumulative_nq_Y)-0.1*cdf_axis_length_Y]

tt_knots_z=[min(cumulative_nq_Z)+0.1*cdf_axis_length_Z,min(cumulative_nq_Z)+0.25*cdf_axis_length_Z,np.mean(cumulative_nq_Z),\
max(cumulative_nq_Z)-0.25*cdf_axis_length_Z,max(cumulative_nq_Z)-0.1*cdf_axis_length_Z]


print 'Knots values for LSQ Spline in CDF = ',tt_knots_z
#t_knots=[(min(x0_Z)+np.mean(x0_Z))/2,np.mean(x0_Z),(max(x0_Z)+np.mean(x0_Z))/2]
#print t_knots
#ff_X = interpolate.LSQUnivariateSpline(cumulative_nq_X, xx_0_X,tt_knots_x)
#ff_Y = interpolate.LSQUnivariateSpline(cumulative_nq_Y, xx_0_Y,tt_knots_y)
#ff_Z = interpolate.LSQUnivariateSpline(cumulative_nq_Z, xx_0_Z,tt_knots_z)


ff_Z = interp1d(cumulative_nq_Z, xx_0_Z)
ff_Y = interp1d(cumulative_nq_Y, xx_0_Y)
ff_X = interp1d(cumulative_nq_X, xx_0_X)


Num_Of_Target_Particles=int(NumberOfParticlesPerSlice*DensityFactor)
 
#density_Z=np.linspace(min(cumulative_nq_Z),max(cumulative_nq_Z),Num_Of_Target_Particles)
density_Z=np.random.uniform(low=min(cumulative_nq_Z), high=max(cumulative_nq_Z), size=(Num_Of_Target_Particles))    
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

#from scipy.interpolate import Rbf

#rbfi_PX = Rbf(mA_X, mA_Y, mA_Z,mA_PX)
#rbfi_PY = Rbf(mA_X, mA_Y, mA_Z,mA_PY)
#rbfi_PZ = Rbf(mA_X, mA_Y, mA_Z,mA_PZ)

#Full_PX = rbfi_PX(Full_X, Full_Y, Full_Z)
#Full_PY = rbfi_PY(Full_X, Full_Y, Full_Z)
#Full_PZ = rbfi_PZ(Full_X, Full_Y, Full_Z)

Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
#FULL_NoParticles=interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_WGHT.ravel(),(Full_X, Full_Y, Full_Z), method='linear')

No_Particles_Per_Record=np.zeros((len(Full_Z),1))

No_Particles_Per_Record[0:len(Full_Z)]=int(round((TotalNumberOfElectrons/len(Full_Z))))


x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,No_Particles_Per_Record.flat]).T
#x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,FULL_NoParticles.flat]).T

#x_px_y_py_z_pz_NE_noNaN=x_px_y_py_z_pz_NE[~np.isnan(x_px_y_py_z_pz_NE).any(axis=1)]

print 'Final charge of particles = ',np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)
#plt.xlabel('Z')
#plt.ylabel('Density')
#plt.plot(x0_Z,y0_Z,'ro')
#plt.ticklabel_format(style='sci',scilimits=(0,0))
#plt.grid(True)
#plt.plot(xnew_Z_plt,ynew_Z_plt)


output_file=tables.open_file(file_name_base+'_FR2.si5','w')

# Create hdf5 file

my_dpi=96
plt.figure(figsize=(1500/my_dpi, 800/my_dpi), dpi=my_dpi)

plt.title('Density profiles')
plt.xlabel('Z')
plt.ylabel('Density')
plt.hist(Full_Z,weights=No_Particles_Per_Record,bins=number_of_bins,label='Density histogram of new data set')
plt.plot(x0_Z,y0_Z,'yo',label='Initial data density values')
plt.ticklabel_format(style='sci')
plt.grid(True)
plt.plot(xnew_Z_plt,ynew_Z_plt,linewidth=4, color='red',label='Fitted density profile')
plt.legend(loc='upper right')

plt.savefig('Final_plot.png',dpi=300)
plt.show()


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

