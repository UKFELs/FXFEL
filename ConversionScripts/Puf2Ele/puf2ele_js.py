# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:52:33 2015

@author: piotrt
"""
# The script is converting Puffin ouput into Elegant format
# Because Elegant requires particles to have same charge
# the script groups particles into fatter particles of same charge.
# It is done in brutal ripping way - see the loop responsible for this.
# Of cours better approach is CRUCIAL to have real results.
# The idea is proposed by Jonathan Smith - just for testing purposes.

import numpy as np
import h5py
import tables
import os
import sys

# Some constants
ERM=0.51099906 #Electron Rest Mass
me=9.11e-31
e_ch=1.602e-19
c=3e+08

# Read the name of the file - step number
if len(sys.argv)==2:
   step_name_in=sys.argv[1]
   print 'Processing file:', step_name_in
else:
   print 'Usage: puf2ele <TimeStep Number e.g. 750> \n'
   sys.exit(1)  

# Assign name of the files basing on given step number
print 'Merging Puffin output files...'
step_no=str(step_name_in)
xdatfile=step_no+'_XDataFile.dat'
ydatfile=step_no+'_YDataFile.dat'
z2datfile=step_no+'_Z2DataFile.dat'
RE_datfile=step_no+'_RE_PPerpDataFile.dat'
IM_datfile=step_no+'_IM_PPerpDataFile.dat'
Gamma_datfile=step_no+'_GammaDataFile.dat'
Chi_datfile='ChiDataFile.dat'
ZDataFile=step_no+'_ZDataFile.dat'
outfile=step_no+'_all.sdds'

final_output=step_no+'_all.sdds'


# Use sddsxref to merge files into one large file
os.system("sddsxref %s %s %s %s %s %s %s %s"\
%(xdatfile,ydatfile,z2datfile,RE_datfile,IM_datfile,Gamma_datfile,Chi_datfile,outfile))

# Convert above SDDS file into HDF to easier use in Python
print 'Converting SDDS Puffin output to HDF5...'

final_output_hdf=step_no+'_all.h5'
os.system("sdds2hdf %s %s"\
%(final_output,final_output_hdf))
os.remove(final_output)

# Load data from above created HDF5 file

print 'Loading data from Puffin HDF5...'

f=h5py.File(final_output_hdf,'r')

X_f=np.array(f['/page1/columns/X'])
Y_f=np.array(f['/page1/columns/Y'])
Z2_f=np.array(f['/page1/columns/Z2'])

Gamma_f=np.array(f['/page1/columns/Gamma'])
IM_PPerp_f=np.array(f['/page1/columns/IM_PPerp'])
RE_PPerp_f=np.array(f['/page1/columns/RE_PPerp'])
s_chi_bar_f=np.array(f['/page1/columns/s_chi_bar'])
f.close()
os.remove(final_output_hdf)

# Another constants
me=9.11E-31
c=3E+08
os.system("sdds2hdf ParamDataFile.dat param.h5")

f2=h5py.File('param.h5','r')
Lc_f=np.array(f2['/page1/parameters/Lc'])
Lg_f=np.array(f2['/page1/parameters/Lg'])
a_u_f=np.array(f2['/page1/parameters/aw'])
gamma_0_f=np.array(f2['/page1/parameters/gamma_r'])
npk_bar=np.array(f2['/page1/parameters/npk_bar'])
f2.close()
os.remove('param.h5')

print 'Read constants...:'
print 'Npk_bar= ',npk_bar[0]
print 'Lc = ',Lc_f[0]
print 'Lg = ',Lg_f[0]
print 'aw = ',a_u_f[0]
print 'gamma_r = ',gamma_0_f[0]

os.system("sdds2hdf %s Z.h5"\
%(ZDataFile))
z_file=h5py.File('Z.h5','r')
z_param_f=np.array(z_file['/page1/columns/Z'])
z_file.close()
os.remove('Z.h5')

t=z_param_f[0]/c
print 'Z = ',z_param_f[0]

z_bar=z_param_f[0]
print 'z_bar = ', z_bar

# Convert data from Puffin to Elegant using equations from Puffin paper.
print 'Processing data...'

z=Lg_f[0]*z_bar

S=(Lc_f[0]*Z2_f[:])+z



X=np.sqrt(Lg_f[0]*Lc_f[0])*X_f[:]
Y=np.sqrt(Lg_f[0]*Lc_f[0])*Y_f[:]

Gamma=Gamma_f[:]*gamma_0_f[0]
P=me*me*c*c((Gamma**2)-1)

Px=RE_PPerp_f[:]*me*c*a_u_f[0]
Py=-1*IM_PPerp_f[:]*me*c*a_u_f[0]
Pz=np.sqrt(P**2-Px**2-Py**2)

Beta_z=Pz/P

xp=Px/Pz
yp=Py/Pz
t=z/(c*Beta_z)


# Create array with all necessary data
full_array=np.column_stack((X,Y,S,Px,Py,Gamma,s_chi_bar_f))
# Create HDF file and SDDS file - final outputs.
FinalHDF_File=step_no+'_full_array.h5'
FinalSDDS_File=step_no+'_full_array.sdds'

print 'Creating final HDF5 output...'
output_file=tables.open_file(FinalHDF_File,'w')

# Create group for particles, don't add metadata yet ! Remove some not necessary metadata
ParticleGroup=output_file.create_array('/','spatialPositions',full_array)
#boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
#boundsGroup._v_attrs.vsType='limits'
#boundsGroup._v_attrs.vsKind='Cartesian'
#timeGroup=output_file.create_group('/','time','time')
#timeGroup._v_attrs.vsType='time'
del ParticleGroup._v_attrs.CLASS
del ParticleGroup._v_attrs.VERSION
del ParticleGroup._v_attrs.FLAVOR
del ParticleGroup._v_attrs.TITLE
#ParticleGroup._v_attrs.vsType='variableWithMesh'
#ParticleGroup._v_attrs.vsTimeGroup='time'
#ParticleGroup._v_attrs.vsNumSpatialDims = 3
#ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
#ParticleGroup._v_attrs.vsLabels='x,y,s,Px,Py,Gamma,Chi'
output_file.close()

# Convert HDF into SDDS
os.system("hdf2sdds %s %s -dataset=spatialPositions"\
%(FinalHDF_File,FinalSDDS_File))

print 'Creating SDDS output...'
# Add proper names of columns in SDDS file
os.system("sddsconvert %s -rename=columns,D1=x,D2=y,D3=s,D4=Px,D5=Py,D6=Gamma,D7=Chi"\
%(FinalSDDS_File))


os.remove(FinalSDDS_File+'~')
# Add metadata to HDF file
print 'Adding some magic...'
output_file=tables.open_file(FinalHDF_File,'a')
ParticleGroup=output_file.root.spatialPositions
output_file=tables.open_file(step_no+'_full_array.h5','a')
boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
boundsGroup._v_attrs.vsType='limits'
boundsGroup._v_attrs.vsKind='Cartesian'
timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'
ParticleGroup._v_attrs.vsType='variableWithMesh'
ParticleGroup._v_attrs.vsTimeGroup='time'
ParticleGroup._v_attrs.vsNumSpatialDims = 3
ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
ParticleGroup._v_attrs.vsLabels='x,y,s,Px,Py,Gamma,Chi'
output_file.close()
print 'Created SDDS: ',FinalSDDS_File
print 'Created HDF5: ',FinalHDF_File
#sys.exit(1)

# This is EXPERIMENTAL
# Intends to group particles in such manner they have same weight
# Uses brutal method. Just counts to desired charge and creates particle...
# Well, it just for testing purposes.

particles=full_array
# Find maximum charge and multiply it by 10.

max_charge=10*np.max(particles[:,6])
print 'Maximum charge= ',max_charge
num_required_particles=10000
total_charge=np.sum(particles[:,6]) # sum charges
# Find total charge
print 'Total charge= ', total_charge
sum_ch=0
newparticles=[] # new array
charges=[]
print newparticles
#charge_required=total_charge/float(num_required_particles) # ensure float division
charge_required=max_charge

#charge_required=10000
print 'Charge required= ', charge_required
print 'Loop over ',np.shape(particles)[0],' number of records...'


# Loop over particles and group them
for i in range(np.shape(particles)[0]-1):
  sum_ch=sum_ch+particles[i,6]

  if (sum_ch)>(charge_required):

     newparticles.append(particles[i])
     charges.append(sum_ch)
#    np.hstack((newparticles, particles[i,0:5]))
     sum_ch=0     
#    charge_required+=total_charge/float(num_required_particles)
     
# Create array with new particles - less than original but with same weight.

newparticles=np.asarray(newparticles)

#print newparticles
part_and_charge=np.column_stack((newparticles,charges))
total_charge_new=np.sum(part_and_charge[:,7])
charge_diff=(total_charge_new/total_charge)*100

# Print the difference in original total charge and current total charge 
print 'Charge factor= ',charge_diff,'%'

# Save the data into HDF file with 'js' at beginning of the name
# No SDDS creation yet <- To be done same way as above.

output_file2=tables.open_file('js'+FinalHDF_File,'w')
#ParticleGroup=output_file2.create_array('/','spatialPositions',particles)
ParticleGroup2=output_file2.create_array('/','GroupedCharge',part_and_charge)
boundsGroup2=output_file2.create_group('/','globalGridGlobalLimits','')
boundsGroup2._v_attrs.vsType='limits'
boundsGroup2._v_attrs.vsKind='Cartesian'
timeGroup2=output_file2.create_group('/','time','time')
timeGroup2._v_attrs.vsType='time'
ParticleGroup2._v_attrs.vsType='variableWithMesh'
ParticleGroup2._v_attrs.vsTimeGroup='time'
ParticleGroup2._v_attrs.vsNumSpatialDims = 3
ParticleGroup2._v_attrs.vsLimits='globalGridGlobalLimits'
ParticleGroup2._v_attrs.vsLabels='x,y,s,Px,Py,Gamma,Chi,Scaled_Charge'
output_file2.close()
