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
# The output file is saved as HDF5 file with VizSchema metadata added

#############################################################


#################################################################################
#################################################################################
#DISABLE WARNINGS - COMMENT IF YOU WISH TO SEE THE WARNINGS FROM BELOW CODE !!!!!
#################################################################################
#################################################################################
import warnings
warnings.filterwarnings('ignore')
#################################################################################
#################################################################################



import numpy as np
import tables
import sys
import time
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import datetime





if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SU2CDF <Input File Name>\n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()

import PARAMS_CDF

#  Set default parameters and read from parameters file

try:
    k_u=PARAMS_CDF.k_u
except:
    k_u = 228.4727
try:
    a_u=PARAMS_CDF.a_u
except:
    a_u = 1.0121809
try:
    DensityFactor=PARAMS_CDF.DensityFactor
except:
    DensityFactor = 10
try:
    SlicesMultiplyFactor=PARAMS_CDF.SlicesMultiplyFactor
except:
    SlicesMultiplyFactor = 20
#try:
#    binnumber=PARAMS_CDF.DensitySampling
#except:
#    binnumber = 15
try:    
    binnumber_X=PARAMS_CDF.X_DensitySampling 
except:
    binnumber_X = 50
try:
    binnumber_Y=PARAMS_CDF.Y_DensitySampling
except:
    binnumber_Y = 50
try:
    binnumber_Z=PARAMS_CDF.Z_DensitySampling
except:
    binnumber_Z = 25
try:
    S_factor=PARAMS_CDF.BeamStretchFactor
except:
    S_factor = 0.0
try:
    NumShapeSlices=PARAMS_CDF.ShapeSampling
except:
    NumShapeSlices = 100

    
if NumShapeSlices < 0:
    NumShapeSlices = 50
print 'user sampling = ',PARAMS_CDF.ShapeSampling


    

#==============================================================================
print 'User defined parameters:'
print 'k_u = ',k_u
print 'a_u = ',a_u
print 'Density Factor = ',DensityFactor
print 'Slices per wavelnegth = ',SlicesMultiplyFactor
#print 'Particle density samples = ',binnumber
print 'Density sampling in X = ',binnumber_X
print 'Density sampling in Y = ',binnumber_Y
print 'Current / Density sampling in Z =' ,binnumber_Z
print 'Stretching factor in Z = ',S_factor
#print 'Shape sampling number = ',NumShapeSlices
#==============================================================================


start = time.time()
def elapsed():
    return time.time() - start


# Pre-allocate array and load data into array



f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()

#Filter the particles with z values



#*************************************************************
# The below section calculates size of the bin according
# to value of Lc (binnumbers=total_length/Lc)
# USER DATA - MODIFY ACCORDING TO REQUIREMENTS
Pi=np.pi                    # Pi number taken from 'numpy' as more precise than just 3.1415
c=299792458.0               # Speed of light
m=9.11e-31                  # mass of electron
e_0=8.854E-12               # vacuum permitivity
e_ch=1.602e-19              # charge of one electron

#*************************************************************
#Filter the particles
#pz_low=615.0
#pz_high=620.0
#z_low=0.0039450
#mA_X=Particles[:,0][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_Y=Particles[:,2][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_Z=Particles[:,4][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_PX=(m*c)*Particles[:,1][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_PY=(m*c)*Particles[:,3][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_PZ=(m*c)*Particles[:,5][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]
#mA_WGHT=Particles[:,6][(Particles[:,5]>=pz_low) & (Particles[:,5]<pz_high) & (Particles[:,4]>z_low)]

mA_X = Particles[:,0]
mA_Y = Particles[:,2]    
mA_Z = Particles[:,4]
# Descale the momentum units from p/mc to SI - to keep compatibility with rest of calculations
mA_PX = Particles[:,1]*(m*c)
mA_PY = Particles[:,3]*(m*c)  
mA_PZ = Particles[:,5]*(m*c)
mA_WGHT = Particles[:,6]

# The below section calculate some initial data - 4*Pi*Rho is the one mose desired
xyz = np.vstack([mA_X,mA_Y,mA_Z]).T
size_x=max(mA_X)-min(mA_X)
size_y=max(mA_Y)-min(mA_Y)
size_z=max(mA_Z)-min(mA_Z)
# This is number of bins just to calculate initial data - don't change if not sure

#cube_volume=(size_x*size_y*size_z)/float(binnumber**3)
#print 'Volume of sample cube = ',cube_volume
print 'Size of sample X,Y,Z = ', size_x,size_y,size_z
#H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=mA_WGHT.flat)
#n_p=float(np.amax(H))/cube_volume
#print 'Particle density N_p = ',n_p


# Calculate some needed values: Omega,Rho,Lc,Lambda_U,Lambda_R,
p_tot=np.sqrt((mA_PX[:]**2)+(mA_PY[:]**2)+(mA_PZ[:]**2))
gamma=(np.sqrt(1+(p_tot/(m*c))**2))
gamma_0=np.mean(gamma)


#omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))
#rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))





lambda_u=(2.0*Pi)/k_u
#print 'Rho = ',rho
# Calculated for planar undulator -> (1+(a^2)/2)
lambda_r=(lambda_u/(2.0*gamma_0**2.0))*(1+(a_u**2.0)/2.0)
#lambda_r=(lambda_u/(2*gamma_0**2))*(1+a_u**2)
#Lc=lambda_r/(4.0*Pi*rho)
print 'lambda_r = ',lambda_r
#print 'L_c = ',Lc
# End of inital data calculations

#*************************************************************

# Stack the particles and charge in common array to allow sorting of data along Z axis
#xyzW = np.vstack((mA_X.flat,mA_Y.flat,mA_Z.flat,mA_WGHT.flat)).T
# Sort the data along Z axis
#xyzW=xyzW[xyzW[:,2].argsort()]


#Calculate total charge
TotalNumberOfElectrons=np.sum(mA_WGHT)

# Assign the sorted data to separate arrays
#m_X=xyzW[:,0].flat
#m_Y=xyzW[:,1].flat
#m_Z=xyzW[:,2].flat
#m_WGHT=xyzW[:,3].flat
NumberOfSourceParticles=len(mA_X)
InitialParticleCharge=TotalNumberOfElectrons*e_ch
# Print some user useful informations

print 'Initial charge of particles = ',InitialParticleCharge
print 'Total number of electrons = ',TotalNumberOfElectrons
print 'Number of source macroparticles = ',len(mA_X)
print 'Electrons/macroparticles in source data = ',round(TotalNumberOfElectrons/len(mA_X))
print 'Desired Electrons/macroparticles in output data = ',int(round(TotalNumberOfElectrons/(len(mA_X)*DensityFactor)))
print 'Macroparticles in job = ',NumberOfSourceParticles


# Merge XZ and YZ into common array to allow creation of 2D histogram data for XZ and YZ planes          
m_Xm_Z=np.vstack((mA_X.flat,mA_Z.flat)).T
m_Ym_Z=np.vstack((mA_Y.flat,mA_Z.flat)).T

# Set the factor to extend histogram with ZERO values to smoothen the edges. Set to 0 if not needed.
# The value of 0.15 means that the histogram will grow 30% in each direction (from -1.30*size to +1.13*size)

#S_factor=0.015

# Create histogram for Z direction and stretch it using S_factor
Hz, edges_Z = np.histogramdd(mA_Z, bins = binnumber_Z,range=((min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=mA_WGHT)
edges_YZ = np.histogramdd(m_Ym_Z, bins = (binnumber_Y,binnumber_Z),normed=False,weights=mA_WGHT)

# Crate histogram for XZ and YZ planes and stretch it using S_factor variable
HxHz,edges_XZ = np.histogramdd(m_Xm_Z, bins = (binnumber_X,binnumber_Z),range=((min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=mA_WGHT)
HyHz,edges_YZ = np.histogramdd(m_Ym_Z, bins = (binnumber_Y,binnumber_Z),range=((min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=mA_WGHT)

# Calculate peak current
#Peak_Current=((np.max(Hz)*e_ch)/(size_z/binnumber_Z))*c
#print 'Peak current = ',Peak_Current

# Peak concetration assumes circular transverse shape and averages over x/y
#Transverse_area=(np.pi*((size_x+size_y)/2.0)**2.0)
#print 'Peak concentration = ',Peak_Current/c/e_ch/(Transverse_area)

# Gaussian smoothing of denisty profiles in X/Z and Y/Z
#Hz=ndimage.gaussian_filter(Hz,1.0)
#HxHz=ndimage.gaussian_filter(HxHz,1.0)
#HyHz=ndimage.gaussian_filter(HyHz,1.0)


# Initiate empty array for data points in histogram (move the histogram to array like: Value,X,Y,Z)
XZarr=np.zeros(((len(edges_XZ[0])-1)*(len(edges_XZ[1])-1),3))
YZarr=np.zeros(((len(edges_YZ[0])-1)*(len(edges_YZ[1])-1),3))


#*******Interpolate density*************

# Z axis interpolation 
# Create equispaced points 
x0_Z = np.linspace(0.5*(edges_Z[0][0]+edges_Z[0][1]),0.5*(edges_Z[0][binnumber_Z]+edges_Z[0][binnumber_Z-1]),binnumber_Z)
y0_Z = Hz



# If user want to use LSQ interpolation then unhash next three lines and comment the line where RBF is used 
#z_hst_lngth=np.max(x0_Z)-np.min(x0_Z)
#t_knots_z=np.linspace(np.min(x0_Z)+0.15*z_hst_lngth,np.max(x0_Z)-0.15*z_hst_lngth,13)
#f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)


#Smoothen the data
#from statsmodels.nonparametric.smoothers_lowess import lowess
#smoothing_factor=0.01
#y0_Z_smth = lowess(y0_Z, x0_Z, frac=smoothing_factor)
#print y0_Z_smth

# Use RBF interpolation for Z-axis, hash next lines and unhash 3 lines for LSQ interpolation above  
#f_Z = interpolate.Rbf(y0_Z_smth[:,0], y0_Z_smth[:,1])
f_Z = interpolate.UnivariateSpline(x0_Z, y0_Z)


#****Below is just for plotting
m_Z_plt=np.linspace(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z,100)
plt.plot(m_Z_plt,f_Z(m_Z_plt))
plt.show()


#****End of plotting Z axis density profile

# Print the value of the bins with non-zero values in histogram
# This tells us how many histogram values have charge and also informs us how many zero-values we added
# with S_factor variable
Non_Zero_Z=float(np.count_nonzero(Hz))
print 'Non-zero histogram values in Z = ',Non_Zero_Z

# Convert XZ/YZ density histograms to XZ_Density/YZ_Density arrays (move the histogram to array like: Value,X,Y,Z)
for zz in range(1,len(edges_XZ[1])):
  for xx in range(1,len(edges_XZ[0])):
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),0]=(edges_XZ[0][xx]+edges_XZ[0][xx-1])*0.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),1]=(edges_XZ[1][zz]+edges_XZ[1][zz-1])*0.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),2]=HxHz[xx-1,zz-1]

for zz in range(1,len(edges_YZ[1])):
  for yy in range(1,len(edges_YZ[0])):
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),0]=(edges_YZ[0][yy]+edges_YZ[0][yy-1])*0.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),1]=(edges_YZ[1][zz]+edges_YZ[1][zz-1])*0.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),2]=HyHz[yy-1,zz-1]

# Calculate the number of slices using the SlicesMultiplyFactor and 4*Pi*rho calculated before
# Do NOT change this line unless you want to set some number of slices not binded with 4*Pi*rho

NumberOfSlices=int(SlicesMultiplyFactor*((max(mA_Z)+S_factor*size_z)-(min(mA_Z)-S_factor*size_z))/(lambda_r))
print 'Number of slices = ',NumberOfSlices

Num_Of_Slice_Particles=int(NumberOfSourceParticles*DensityFactor/NumberOfSlices)
print 'Number of particles in each slice = ',Num_Of_Slice_Particles

# Calculate the step size - this is just size of samples divided over number of calculated Slices
Step_Size=((max(mA_Z)+S_factor*size_z)-(min(mA_Z)-S_factor*size_z))/NumberOfSlices


#*** INTERPOLATE XZ AND YZ PLANES USING 2D FUNCTION
# Calculate the length of X,Y,Z histogram for fitting
x_hst_lngth=np.max(XZarr[:,0])-np.min(XZarr[:,0])
y_hst_lngth=np.max(YZarr[:,0])-np.min(YZarr[:,0])
z_hst_lngth=np.max(XZarr[:,1])-np.min(XZarr[:,1])

# Calculate knots (t) needed for LSQBivariateSpline
t_XZ=np.linspace(np.min(XZarr[:,0])+0.1*x_hst_lngth,np.max(XZarr[:,0])-0.1*x_hst_lngth,15)
t_YZ=np.linspace(np.min(YZarr[:,0])+0.1*y_hst_lngth,np.max(YZarr[:,0])-0.1*y_hst_lngth,15)
t_ZZ=np.linspace(np.min(XZarr[:,1])+0.1*z_hst_lngth,np.max(XZarr[:,1])-0.1*z_hst_lngth,15)


# Interpolate using LSQBivariateSpline, hash if want to use Interp2D
f_Dens_XZ=interpolate.LSQBivariateSpline(XZarr[:,0].ravel(), XZarr[:,1].ravel(),XZarr[:,2].ravel(),t_XZ,t_ZZ)
f_Dens_YZ=interpolate.LSQBivariateSpline(YZarr[:,0].ravel(), YZarr[:,1].ravel(),YZarr[:,2].ravel(),t_YZ,t_ZZ)

# Interpolate using Interp2D - unhash if want to use instead of LSQ
#f_Dens_XZ=interpolate.interp2d(XZarr[:,0].ravel(), XZarr[:,1].ravel(),XZarr[:,2].ravel())
#f_Dens_YZ=interpolate.interp2d(YZarr[:,0].ravel(), YZarr[:,1].ravel(),YZarr[:,2].ravel())

#*****Plot the density profile XZ,YZ
# Create data for plot

PLT_X=np.linspace(min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x,200)
PLT_Y=np.linspace(min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y,200)
PLT_Z=np.linspace(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z,200)
plt.pcolormesh(PLT_Z, PLT_X,f_Dens_XZ(PLT_X, PLT_Z))
plt.show()
plt.pcolormesh(PLT_Z, PLT_Y,f_Dens_YZ(PLT_Y, PLT_Z))
plt.show()

#*****End of plot

# Initiate data for fitting density profile in each loop
Fit_Particle_Number=1000

New_X=np.linspace(min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x,Fit_Particle_Number)
New_Y=np.linspace(min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y,Fit_Particle_Number)
New_Z=np.zeros(Fit_Particle_Number)

# Initiate empty array for Z positions or particles
density_Z=np.zeros(Num_Of_Slice_Particles)

#Initiate random placed particles with values 0-1 which will be next projected using CDF onto X,Y positions
# The 0-1 is CDF value not the X,Y


from math import log, floor, ceil, fmod
import numpy as np

def halton(dim, nbpts):
    h = np.empty(nbpts * dim)
    h.fill(np.nan)
    p = np.empty(nbpts)
    p.fill(np.nan)
    P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    lognbpts = log(nbpts + 1)
    for i in range(dim):
        b = P[i]
        n = int(ceil(lognbpts / log(b)))
        for t in range(n):
            p[t] = pow(b, -(t + 1) )

        for j in range(nbpts):
            d = j + 1
            sum_ = fmod(d, b) * p[0]
            for t in range(1, n):
                d = floor(d / b)
                sum_ += fmod(d, b) * p[t]

            h[j*dim + i] = sum_

    return h.reshape(nbpts, dim)

#density_Y=np.random.uniform(low=0, high=1, size=(Num_Of_Slice_Particles))
#density_X=np.random.uniform(low=0, high=1, size=(Num_Of_Slice_Particles))
HALTON_XY=halton(2,Num_Of_Slice_Particles)
#print 'Halton shape is = ',np.shape(HALTON_XY)
density_X=HALTON_XY[:,0]
density_Y=HALTON_XY[:,1]

Slice_Ne=np.zeros(Num_Of_Slice_Particles)

# Calculate the min/max values for x/y along z-axis (outer shape)

minz=np.min(mA_Z)
maxz=np.max(mA_Z)
step=(maxz-minz)/NumShapeSlices

# Generate step for rescaled lenght (s_factor)
step_s_factor=((maxz+(S_factor*size_z))-(minz-(S_factor*size_z)))/NumShapeSlices
mmax_X=np.zeros(0)
mmin_X=np.zeros(0)
mmax_Y=np.zeros(0)
mmin_Y=np.zeros(0)
mm_Z=np.zeros(0)



#def Weighted_Stats(values, weights):
#    Return the weighted average and standard deviation

#    average = np.average(values, weights=weights)
#    variance = np.average((values-average)**2.0, weights=weights)  # Fast and numerically precise
#    return (np.sqrt(variance))

# Create interpolated function which describes outer boundaries of initial electron beam
for i in range(0,NumShapeSlices):
#    print minz+step*(i),' ',i
#   Check if there are particles in shape slice
    ShapeParticles = mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]
#    print len(ShapeParticles)
    if len(ShapeParticles)>3:
#        mean_x=np.average(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])
#        std_x=Weighted_Stats(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))],weights=mA_WGHT[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])
#        std_x=np.std(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])         
#        mmax_X=np.append(mmax_X,mean_x+std_x)
#        mmin_X=np.append(mmin_X,mean_x-std_x)
    
#        mean_y=np.average(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])
#        std_y=Weighted_Stats(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))],weights=mA_WGHT[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])
#        std_y=np.std(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))])         
#        mmax_Y=np.append(mmax_Y,mean_y+std_y)
#        mmin_Y=np.append(mmin_Y,mean_y-std_y)
        mmax_X=np.append(mmax_X,np.max(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmin_X=np.append(mmin_X,np.min(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmax_Y=np.append(mmax_Y,np.max(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmin_Y=np.append(mmin_Y,np.min(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
    #    mm_Z[i]=0.5*((minz+step*(i))+(minz+step*(i+1)))
    #   Project the shape of beam on rescaled length (s_factor) 
        mm_Z=np.append(mm_Z,0.5*(((minz-(S_factor*size_z))+step_s_factor*(i))+((minz-(S_factor*size_z))+step_s_factor*(i+1))))

#f_mmax_X=interpolate.UnivariateSpline(mm_Z,mmax_X)
#f_mmin_X=interpolate.UnivariateSpline(mm_Z,mmin_X)
#f_mmax_Y=interpolate.UnivariateSpline(mm_Z,mmax_Y)
#f_mmin_Y=interpolate.UnivariateSpline(mm_Z,mmin_Y)

# Choose below between cubic or linear approximation of shape
#Smoothen the data
#from statsmodels.nonparametric.smoothers_lowess import lowess
#smoothing_factor=0.05
#mmax_X_smth = lowess(mmax_X, mm_Z, frac=smoothing_factor)
#mmin_X_smth = lowess(mmin_X, mm_Z, frac=smoothing_factor)
#mmax_Y_smth = lowess(mmax_Y, mm_Z, frac=smoothing_factor)
#mmin_Y_smth = lowess(mmin_Y, mm_Z, frac=smoothing_factor)
#==============================================================================
f_mmax_X=interpolate.Rbf(mm_Z,mmax_X,function='cubic')
f_mmin_X=interpolate.Rbf(mm_Z,mmin_X,function='cubic')
f_mmax_Y=interpolate.Rbf(mm_Z,mmax_Y,function='cubic')
f_mmin_Y=interpolate.Rbf(mm_Z,mmin_Y,function='cubic')
#==============================================================================
   
#*** Procedure for placing electrons in each slice according to calculated CDF

# Create values for maximum values of X/Y of particle set
OldRange_X=max(New_X)-min(New_X)
OldRange_Y=max(New_Y)-min(New_Y)
OldMin_X=min(New_X)
#OldMax_X=max(New_X)
OldMin_Y=min(New_Y)
#OldMax_Y=max(New_Y)


SliceNumber=np.zeros(Num_Of_Slice_Particles)

def SliceCalculate(slice_number,NumShapeSlicesIN):

# Calculate value (in meters) of current slice
    Z_Slice_Value=(slice_number*Step_Size)+(np.min(mA_Z)-S_factor*size_z)
    density_Z[:]=Z_Slice_Value
    SliceNumber[:]=slice_number

# Scale the range of new particles to new particles set (keep the new particles within original shape of beam)
    if NumShapeSlicesIN > 0:
        NewRange_X=f_mmax_X(Z_Slice_Value)-f_mmin_X(Z_Slice_Value)
        NewRange_Y=f_mmax_Y(Z_Slice_Value)-f_mmin_Y(Z_Slice_Value)
        New_Xl=(((New_X-OldMin_X)*NewRange_X)/OldRange_X)+f_mmin_X(Z_Slice_Value)
        New_Yl=(((New_Y-OldMin_Y)*NewRange_Y)/OldRange_Y)+f_mmin_Y(Z_Slice_Value)
    else:
        New_Xl=New_X
        New_Yl=New_Y
# Interpolate density curve for current slice and remove values below 0
      
    Dens_XZ=f_Dens_XZ(New_X,Z_Slice_Value)
    Dens_YZ=f_Dens_YZ(New_Y,Z_Slice_Value)
    Dens_XZ=Dens_XZ.clip(min=0)
    Dens_YZ=Dens_YZ.clip(min=0)   

# Check if the density is above 0, we don't want to create slices
# with electrons charge = 0
# Note that due to used algorithm macroparticle might have charge less than
# charge of single electron - Puffin accepts this type of data
    if (np.sum(Dens_XZ) > 0 and np.sum(Dens_YZ) > 0): 

# Scale density profile so that the max value is 1
        Dens_XZ=Dens_XZ/np.sum(Dens_XZ)
        Dens_YZ=Dens_YZ/np.sum(Dens_YZ)
   
# Calculate CDF
        cumulative_XZ=np.cumsum(Dens_XZ)
        cumulative_YZ=np.cumsum(Dens_YZ)
    
# Sort the CDF and remove duplicates [MOT VALID ANYMORE]

#        cumulative_nq_XZ=sorted(set(cumulative_XZ))
#        cumulative_nq_YZ=sorted(set(cumulative_YZ))  
# Create the CDF function
        cumulative_nq_XZ=cumulative_XZ
        cumulative_nq_YZ=cumulative_YZ  


# Create linear array for interpolation of CDF onto new X,Y,Z
        xx_0_XZ=np.linspace(np.min(New_Xl),np.max(New_Xl),len(cumulative_nq_XZ))
        xx_0_YZ=np.linspace(np.min(New_Yl),np.max(New_Yl),len(cumulative_nq_YZ))
# Do NOT use the outer boundary checks !
#        xx_0_XZ=np.linspace(np.min(New_X),np.max(New_X),len(cumulative_nq_XZ))
#        xx_0_YZ=np.linspace(np.min(New_Y),np.max(New_Y),len(cumulative_nq_YZ))

# Create CDF interpolation function using  interp1d 
        ff_XZ = interpolate.interp1d(cumulative_nq_XZ, xx_0_XZ)
        ff_YZ = interpolate.interp1d(cumulative_nq_YZ, xx_0_YZ)
        
# Calculate the charge for current slice taking into account that there were some
# slices with charge 0 added by using S_factor        
        Slice_Ne[:]=Non_Zero_Z*f_Z(Z_Slice_Value)/float(NumberOfSlices)

# If the charge of slice is > 0 then apply calculated above CDF and position the electrons
# according to CDF      
        if (np.sum(Slice_Ne))>0:        
            Full_Xl=ff_XZ(density_X)
            Full_Yl=ff_YZ(density_Y)
            Full_Zl=density_Z       
            Full_Nel=Slice_Ne/len(Slice_Ne)
# Return the values to main program
#        print 'Mean Loop Z =',np.mean(Full_Zl),'Slice No =',slice_number            
    return Full_Xl,Full_Yl,Full_Zl,Full_Nel,SliceNumber

#***End of procedure

#***Start multiprocessing of 'SliceCalculate' procedure
#==============================================================================

def main():
    import multiprocessing
    pool = multiprocessing.Pool()
    if __name__ == '__main__':
        result = []
        for slice_number in range(0,NumberOfSlices):
            ProgressValue=100*float(slice_number)/float(NumberOfSlices)
            print 'Completed = ',ProgressValue,' [%]\r',
            pool.apply_async(SliceCalculate, (slice_number,PARAMS_CDF.ShapeSampling,), callback=result.append)
        pool.close()
        pool.join()
    return result
result=main()

#==============================================================================
#***End of multiprocessing

#Rearrange data from multiprocessing into more convenient arrays

Total_Number_Of_Particles=Num_Of_Slice_Particles*len(result)
print 'Total particles in set = ',Total_Number_Of_Particles

#Inititate empty arrays for positions, momentum and charge of macroparticles 
Full_X=np.zeros(Total_Number_Of_Particles)
Full_PX=np.zeros(0)
Full_Y=np.zeros(Total_Number_Of_Particles)
Full_PY=np.zeros(0)
Full_Z=np.zeros(Total_Number_Of_Particles)
Full_PZ=np.zeros(0)
Full_Ne=np.zeros(Total_Number_Of_Particles)

# Loopt that rearranges array 'results' into separate X,Y,Z,Ne arrays
print 'Starting rearranging array...'
counter=0
for j in range(0,Num_Of_Slice_Particles):
    for i in range(0,len(result)):
        Full_X[counter]=result[i][0][j]
        Full_Y[counter]=result[i][1][j]
        Full_Z[counter]=result[i][2][j]
        Full_Ne[counter]=result[i][3][j]
        counter=counter+1


# Add noise
print 'Adding noise...'
Rand_Z=(Step_Size*(np.random.random(len(Full_Z)) - 0.50))/np.sqrt(Full_Ne)
Full_Z=Full_Z+Rand_Z

# print 'Distance between slices = ',Step_Size
#print 'Average noise shift = ',np.mean((Rand_Z/np.sqrt(Full_Ne)))
#print 'Min and Max noise shift = ',np.min((Rand_Z/np.sqrt(Full_Ne))),' ',np.max((Rand_Z/np.sqrt(Full_Ne)))
#print 'Noise standard deviation = ',np.std((Rand_Z/np.sqrt(Full_Ne)))   
# Interpolate momentum onto new macroparticles using the momentum map from initial data   
print 'Starting to interpolate momentum data... - takes time' 

def Calculate_PX(): 
    Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest',rescale=True)
    return Full_PX
def Calculate_PY():
    Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest',rescale=True)
    return Full_PY
def Calculate_PZ(): 
    Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='linear',rescale=True)
    return Full_PZ
    
# Start three simultaneous processes for griddata interpolation
from multiprocessing.pool import ThreadPool
pool = ThreadPool(processes=3)
async_result_PX = pool.apply_async(Calculate_PX, ()) 
async_result_PY = pool.apply_async(Calculate_PY, ())
async_result_PZ = pool.apply_async(Calculate_PZ, ())
Full_PX = async_result_PX.get()
Full_PY = async_result_PY.get()
Full_PZ = async_result_PZ.get()
# End of interpolation



# Open output file 
output_file=tables.open_file(file_name_base+'_CDF.h5','w')

# Scale units from SI to p = p/mc

Full_PX=Full_PX/(m*c)
Full_PY=Full_PY/(m*c)
Full_PZ=Full_PZ/(m*c)


# Merge all data into one array
x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,Full_Ne.flat]).T

# print 'Add Poisson noise to particle weights...'  (AJTC)
x_px_y_py_z_pz_NE[:,6]=np.random.poisson(x_px_y_py_z_pz_NE[:,6])

# Remove all particles with weights <= zero
x_px_y_py_z_pz_NE=x_px_y_py_z_pz_NE[x_px_y_py_z_pz_NE[:,6] > 0]

# Check for NaN calues in data - happens when using 'linear' option in momentum interpolations (griddata parameter)
NaN_Mask=~np.any(np.isnan(x_px_y_py_z_pz_NE), axis=1)
x_px_y_py_z_pz_NE=x_px_y_py_z_pz_NE[NaN_Mask]



# Rescale the charge of new particle set (needed due to S_factor usage)
ChargeFactor=InitialParticleCharge/np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)
print 'Charge scaling factor = ',ChargeFactor
x_px_y_py_z_pz_NE[:,6]=x_px_y_py_z_pz_NE[:,6]*ChargeFactor

print 'Final charge of particles = ',np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)
print 'Saving the output to files...'

# Create Group in HDF5 and add metadata for Visit
ParticleGroup=output_file.create_array('/','Particles',x_px_y_py_z_pz_NE)
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
now = datetime.datetime.now()
ParticleGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
ParticleGroup._v_attrs.FXFELDensityFactor=DensityFactor
#Close the file
output_file.close()

