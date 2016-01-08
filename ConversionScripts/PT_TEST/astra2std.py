# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:55:37 2015

@author: piotrt
"""
import numpy as np
#set some variables, material constants
import sys

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: astra2std <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()

em=0.51099906e6
me=1.0/em

# Open and load file, then read as lines into array
f2=open(file_name_in,'r')
lines=f2.readlines()
f2.close()
data=[]
for line in lines:
		p=line.split()
		data.append(p)
data=[line.split() for line in lines]

# define array as float
data2=np.array(data)
data2=data2.astype(float)
print(data2.shape)

# assign chosen rows from loaded array
x=data2[:,0]
y=data2[:,1]
z=data2[:,2]
px=data2[:,3]
py=data2[:,4]
pz=data2[:,5]
dumb_zero=data2[:,6]
macrocharge=data2[:,7]

# prepare array - needed for Python to work properly
norm_px=px[:]
norm_py=py[:]
norm_pz=pz[:]


# scale first-reference particle
norm_px[0]=px[0]*me
norm_py[0]=py[0]*me
norm_pz[0]=pz[0]*me 


#pz needs to be recalculated again as function of reference particle
# scale rest of particles momenta using the reference particle
for i in range(1,len(x)):
    z[i]=z[i]+z[0]
    norm_pz[i]=pz[i]*me+pz[0]
    norm_py[i]=py[i]*me
    norm_px[i]=px[i]*me
    
 
# Calculate average values - centre of mass
xavg=np.mean(x)
yavg=np.mean(y)
zavg=np.mean(z)
pxavg=np.mean(norm_px)
pyavg=np.mean(norm_py)
pzavg=np.mean(norm_pz)

#print xavg,yavg,zavg,pxavg,pyavg,pzavg

gam=norm_pz[0]+1

# scale positions using centre of mass
norm_x=(x[:]-xavg)
norm_y=(y[:]-yavg)
norm_z=(z[:]-zavg)


# scale momenta using gamma
norm_px=(norm_px[:]-pxavg)/gam
norm_py=(norm_py[:]-pyavg)/gam
norm_pz=(norm_pz[:]-pzavg)/gam
		
  
output_file=open(file_name_base+'ASTRA.txt','w')

for i in range(len(x)):
#		output_file.write("% .4e" %(norm_x[i]) + \
#                    		" % .4e" %(macrocharge[i]) + "\n")
#for i in range(len(x)):
		output_file.write("% .4E" %(norm_x[i]) + \
		" % .4E" %(norm_y[i]) + \
		" % .4E" %(norm_z[i]) + \
		" % .4E" %(norm_px[i]) + \
		" % .4E" %(norm_py[i]) + \
		" % .4E" %(norm_pz[i]) + \
           " % .4E" %(dumb_zero[i]) + \
           " % .4E" %(macrocharge[i]) + "\n")

output_file.close()
from mpl_toolkits.mplot3d import Axes3D
mr_X=norm_x[0::5]
mr_Y=norm_y[0::5]
mr_Z=norm_z[0::5]
mr_PX=norm_px[0::5]
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(mr_X, mr_Y, mr_Z,c=mr_PX)
#ax.scatter(x,y,z,c=pz)
plt.show()