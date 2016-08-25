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
from scipy import interpolate

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
number_of_bins=100
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
DensityFactor=10
print 'Initial charge of particles = ',TotalNumberOfElectrons*e_ch
print 'Total number of electrons: ',TotalNumberOfElectrons
print 'Number of source macroparticles: ',len(mB_X)
print 'Electrons/macroparticles in source data:',round(TotalNumberOfElectrons/len(mB_X))
print 'Desired Electrons/macroparticles in output data:',int(round(TotalNumberOfElectrons/(len(mB_X)*DensityFactor)))

# Set desired particle number per slice - best to set the number that will give 
# an integer when total number of source particles will be divided over it
# Otherwise your particles will be cut at the end
# Very low number can give weird results

NumberOfSourceParticles=len(mB_X)


print 'Macroparticles in job= ',NumberOfSourceParticles



 
m_X=mB_X[:]
m_Y=mB_Y[:]
m_Z=mB_Z[:]
m_WGHT=mB_WGHT[:]
          
m_Xm_Z=np.vstack((mA_X.flat,mA_Z.flat)).T
m_Ym_Z=np.vstack((mA_Y.flat,mA_Z.flat)).T
print np.shape(m_Xm_Z),np.shape(m_Ym_Z)

binnumber_Z=number_of_bins   
#binnumber_Z=100

# Rescale X/Y bin number to 1/10 of the Z bin number - this allows more smoother data in case of sparse particle mesh
binnumber_X=int(number_of_bins/10)
binnumber_Y=int(number_of_bins/10)

print'Binnumber X,Y,Z = ',binnumber_X,binnumber_Y,binnumber_Z
Hz, edges_Z = np.histogramdd(m_Z, bins = binnumber_Z,normed=False,weights=m_WGHT.flat)
#Hy, edges_Y = np.histogramdd(m_Y, bins = binnumber_Y,normed=False,weights=m_WGHT.flat)
#Hx, edges_X = np.histogramdd(m_X, bins = binnumber_X,normed=False,weights=m_WGHT.flat)


HxHz,edges_XZ = np.histogramdd(m_Xm_Z, bins = (binnumber_X,binnumber_Z),normed=False,weights=m_WGHT.flat)
HyHz,edges_YZ = np.histogramdd(m_Ym_Z, bins = (binnumber_Y,binnumber_Z),normed=False,weights=m_WGHT.flat)


    #import matplotlib.pyplot as plt
    #from scipy.stats import kde
    #data= np.vstack((m_X.ravel(),m_Z.ravel())).T
    #x, y = data.T
    #nbins=20
    #k = kde.gaussian_kde(data.T,weights=mB_WGHT[:])
    #xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    #zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
    #plt.show()

XZarr=np.zeros(((len(edges_XZ[0])-1)*(len(edges_XZ[1])-1),3))
YZarr=np.zeros(((len(edges_YZ[0])-1)*(len(edges_YZ[1])-1),3))


# Interpolate density along Z-axis
#x0_Z = np.linspace(np.min(m_Z),np.max(m_Z),binnumber_Z)
x0_Z = np.linspace(0.5*(edges_Z[0][0]+edges_Z[0][1]),0.5*(edges_Z[0][binnumber_Z]+edges_Z[0][binnumber_Z-1]),binnumber_Z)

y0_Z = Hz

z_hstgrm_length=0.5*(edges_Z[0][binnumber_Z]+edges_Z[0][binnumber_Z-1])-0.5*(edges_Z[0][0]+edges_Z[0][1])

#t_knots_z=[min(x0_Z)+0.1*z_hstgrm_length,min(x0_Z)+0.25*z_hstgrm_length,np.mean(x0_Z),max(x0_Z)-0.25*z_hstgrm_length,max(x0_Z)-0.1*z_hstgrm_length]
t_knots_z=[edges_Z[0][0]+0.25*z_hstgrm_length,(edges_Z[0][binnumber_Z]+edges_Z[0][0])*0.5,edges_Z[0][binnumber_Z]-0.25*z_hstgrm_length]
#t_knots_z=[edges_Z[0][0]+0.25*z_hstgrm_length,edges_Z[0][binnumber_Z]-0.25*z_hstgrm_length]
#f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)
f_Z=interpolate.UnivariateSpline(x0_Z, y0_Z,ext=3)


#f_Z = interpolate.UnivariateSpline(x0_Z, y0_Z,k=5)
#print 0.5*(edges_Z[0][binnumber_Z]+edges_Z[0][binnumber_Z-1]),0.5*(edges_Z[0][0]+edges_Z[0][1])
#print np.min(m_Z),np.max(mA_Z)
plt.plot(x0_Z,f_Z(x0_Z))
plt.show()
# Convert XZ/YZ density histograms to XZ_Density/YZ_Density arrays
for zz in range(1,len(edges_XZ[1])):
  for xx in range(1,len(edges_XZ[0])):
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),0]=(edges_XZ[0][xx]+edges_XZ[0][xx-1])*.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),1]=(edges_XZ[1][zz]+edges_XZ[1][zz-1])*.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),2]=HxHz[xx-1,zz-1]

for zz in range(1,len(edges_YZ[1])):
  for yy in range(1,len(edges_YZ[0])):
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),0]=(edges_YZ[0][yy]+edges_YZ[0][yy-1])*.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),1]=(edges_YZ[1][zz]+edges_YZ[1][zz-1])*.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),2]=HyHz[yy-1,zz-1]

print np.shape(XZarr),np.shape(YZarr)

# Generate random initial slice of X/Y particles set for CDF function - full range (0-1)
# Slices multiply factor is the number how many layers should be within one Lc distance
SlicesMultiplyFactor=10

NumberOfSlices=(number_of_bins+1)*SlicesMultiplyFactor
Num_Of_Slice_Particles=int(NumberOfSourceParticles*DensityFactor/NumberOfSlices)


Step_Size=(np.max(m_Z)-np.min(m_Z))/NumberOfSlices

density_Y=np.random.uniform(low=0, high=1, size=(Num_Of_Slice_Particles))
density_X=np.random.uniform(low=0, high=1, size=(Num_Of_Slice_Particles))

Full_X=np.zeros(0)
Full_PX=np.zeros(0)
Full_Y=np.zeros(0)
Full_PY=np.zeros(0)
Full_Z=np.zeros(0)
Full_PZ=np.zeros(0)
Full_Ne=np.zeros(0)



slice_counter=0
for slice_number in range(0,NumberOfSlices):
       
    Z_Slice_Value=(slice_number*Step_Size)+np.min(m_Z)
    density_Z=np.zeros(Num_Of_Slice_Particles)
    density_Z[:]=(slice_number*Step_Size)+np.min(m_Z)
    ProgressValue=100*((Z_Slice_Value-min(m_Z))/size_z)
#    print 'Z slice value = ',Z_Slice_Value,' m\r',
    print 'Completed = ',ProgressValue,' [%] ',Z_Slice_Value,' [m] \r',

# Interpolate curve density for selected slice


    New_X=np.linspace(min(m_X),max(m_X),100)
    New_Y=np.linspace(min(m_Y),max(m_Y),100)
   
    New_Z=np.zeros(100)
    New_Z=Z_Slice_Value

    XZarr_SL=XZarr[((XZarr[:,1]<=(Z_Slice_Value+(Step_Size*20))) & (XZarr[:,1]>=Z_Slice_Value-(Step_Size*20)))]
    YZarr_SL=YZarr[((YZarr[:,1]<=(Z_Slice_Value+(Step_Size*20))) & (YZarr[:,1]>=Z_Slice_Value-(Step_Size*20)))]


#    Dens_XZ = interpolate.griddata((XZarr_SL[:,0].ravel(), XZarr_SL[:,1].ravel()),XZarr_SL[:,2].ravel(),(New_X, New_Z), method='cubic',fill_value=0)
#    Dens_YZ = interpolate.griddata((YZarr_SL[:,0].ravel(), YZarr_SL[:,1].ravel()),YZarr_SL[:,2].ravel(),(New_Y, New_Z), method='cubic',fill_value=0)

    Dens_XZ = interpolate.griddata((XZarr[:,0].ravel(), XZarr[:,1].ravel()),XZarr[:,2].ravel(),(New_X, New_Z), method='cubic',fill_value=0)
    Dens_YZ = interpolate.griddata((YZarr[:,0].ravel(), YZarr[:,1].ravel()),YZarr[:,2].ravel(),(New_Y, New_Z), method='cubic',fill_value=0)
    
    Dens_XZ=Dens_XZ.clip(min=0)
    Dens_YZ=Dens_YZ.clip(min=0)   

    if (np.sum(Dens_XZ) > 0 and np.sum(Dens_YZ) > 0):     
        Dens_XZ=Dens_XZ/np.sum(Dens_XZ)
        Dens_YZ=Dens_YZ/np.sum(Dens_YZ)
    
  
          
        cumulative_XZ=np.cumsum(Dens_XZ)
        cumulative_YZ=np.cumsum(Dens_YZ)
    
        cumulative_nq_XZ=sorted(set(cumulative_XZ))
        cumulative_nq_YZ=sorted(set(cumulative_YZ))  
    
        xx_0_XZ=np.linspace(np.min(m_X),np.max(m_X),len(cumulative_nq_XZ))
        xx_0_YZ=np.linspace(np.min(m_Y),np.max(m_Y),len(cumulative_nq_YZ))
    
        ff_XZ = interpolate.UnivariateSpline(cumulative_nq_XZ, xx_0_XZ,ext=1)
        ff_YZ = interpolate.UnivariateSpline(cumulative_nq_YZ, xx_0_YZ,ext=1)
        
        Slice_Ne=np.zeros(Num_Of_Slice_Particles)
        
        
   # Charge in slice
        slice_counter=slice_counter+1
        Slice_Ne[:]=(f_Z(slice_counter*Step_Size+np.min(m_Z)))/SlicesMultiplyFactor        
        if (np.sum(Slice_Ne))>=0:        
        #Slice_Ne[:]=(f_Z(slice_number*Step_Size+np.min(m_Z)))/SlicesMultiplyFactor
#           Slice_Ne[:]=1.0
            Full_X=np.append(ff_XZ(density_X),Full_X)
            Full_Y=np.append(ff_YZ(density_Y),Full_Y)
            Full_Z=np.append(density_Z,Full_Z)        
            Full_Ne=np.append(Slice_Ne/len(Slice_Ne),Full_Ne)
    
print ' \r'    
print np.shape(Full_X),np.shape(Full_Y),np.shape(Full_Z),np.shape(Full_Ne)    
 
Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')
Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='nearest')

Rand_Z=Step_Size*(np.random.random(len(Full_Z)) - 0.5)
#print Rand_Z
Full_Z=Full_Z+(Rand_Z/np.sqrt(Full_Ne))

 
output_file=tables.open_file(file_name_base+'_NS.si5','w')

x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,Full_Ne.flat]).T
print 'Final charge of particles = ',np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)  


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



#plt.title('Density profiles')
#plt.xlabel('Z')
#plt.ylabel('Density')
#plt.hist(Full_Z,weights=Full_Ne,bins=number_of_bins,label='Density histogram of new data set')
#plt.grid(True)
#plt.legend(loc='upper right')
#plt.show()

