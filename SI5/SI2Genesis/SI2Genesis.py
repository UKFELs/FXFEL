# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:56:42 2015

@author: piotrt
"""
import tables
import sys
import numpy as np
m=9.11e-31
e_ch=1.602e-19
c=3.0e+8

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SI2Genesis <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()




f=tables.open_file(file_name_in,'r')
Electrons=f.root.Particles.read()

# Assign number of records (macroparticles) as n
n=len(Electrons)

X = Electrons[:,0]
Px = Electrons[:,1]/(m*c)
Y = Electrons[:,2]
Py = Electrons[:,3]/(m*c)
Z = Electrons[:,4]
P = (np.sqrt(Electrons[:,1]**2+Electrons[:,3]**2+Electrons[:,5]**2))/(m*c)
# Write the whole array in MASP input format X,Px,Y,Py,Z,Pz,Weight_of_Particle
# This is done in loop - reading each item from Electrons array

out=open(file_name_base+'_SI2Gen.txt','w')


out.write('# Input Distribution from '+str(file_name_in)+"\n")
out.write('? VERSION = 1.0'+"\n")
out.write('? CHARGE = '+str(np.sum(Electrons[:,6]*e_ch))+"\n")
out.write('? SIZE = '+str(n)+"\n")
out.write('? COLUMNS X PX Y PY Z P'+"\n")


for i in range(n):                                  
                out.write("%.8e" %(X[i]) + \
                " %.8e" %(Px[i]) + \
                " %.8e" %(Y[i]) + \
                " %.8e" %(Py[i]) + \
                " %.8e" %(Z[i]) + \
                " %.8e" %(P[i]) + "\n")                        
                
#Close the files                
out.close()
f.close()
# End of script