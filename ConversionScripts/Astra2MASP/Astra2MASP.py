# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:55:37 2015

@author: piotrt
"""
import numpy as np
import sys
#set some variables, material constants


e_ch=1.602e-19

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: astra2std <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()



# Open and load file, then read as lines into array
file_in=open(file_name_in,'r')
lines=file_in.readlines()
file_in.close()
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
norm_x=x[:]
norm_y=y[:]
norm_z=z[:]

# scale first-reference particle
for i in range(1,len(x)):
#    norm_px[i]=(px[i]+px[0])*me
#    norm_py[i]=(py[i]+py[0])*me
#    norm_pz[i]=(pz[i]+pz[0])*me 
#    norm_x[i]=(x[i]+x[0])
#    norm_y[i]=(y[i]+y[0])
#    norm_z[i]=(z[i]+z[0])
    norm_px[i]=(px[i])*5.36E-28
    norm_py[i]=(py[i])*5.36E-28
    norm_pz[i]=(pz[i]+pz[0])*5.36E-28 
    norm_x[i]=(x[i])
    norm_y[i]=(y[i])
    norm_z[i]=(z[i]+z[0])


norm_px[0]=(px[0])*5.36E-28
norm_py[0]=(py[0])*5.36E-28
norm_pz[0]=(pz[0])*5.36E-28 

# scale positions using centre of mass


	
  
output_file=open(file_name_base+'STD','w')

for i in range(len(x)):
#		output_file.write("% .4e" %(norm_x[i]) + \
#                    		" % .4e" %(macrocharge[i]) + "\n")
#for i in range(len(x)):
		output_file.write("% .4E" %(norm_x[i]) + \
		" % .4E" %(norm_px[i]) + \
		" % .4E" %(norm_y[i]) + \
		" % .4E" %(norm_py[i]) + \
		" % .4E" %(norm_z[i]) + \
		" % .4E" %(norm_pz[i]) + \
#          " % .4E" %(dumb_zero[i]) + \
# Scale nC to C ASTRA uses nC  
           " % .4E" %((macrocharge[i]*(1.E-9))/e_ch) + "\n")

output_file.close()
#from mpl_toolkits.mplot3d import Axes3D
#mr_X=norm_x[0::5]
#mr_Y=norm_y[0::5]
#mr_Z=norm_z[0::5]
#mr_PX=norm_px[0::5]
#import matplotlib.pyplot as plt
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(mr_X, mr_Y, mr_Z,c=mr_PX)
#ax.scatter(x,y,z,c=pz)
#plt.show()