# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:56:42 2015

@author: piotrt
"""
import tables
import sys

m=9.11e-31

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SI_52MASP <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()




f=tables.open_file(file_name_in,'r')
Electrons=f.root.Particles.read()

# Assign number of records (macroparticles) as n
n=len(Electrons)

# Create output file named output.txt
out=open(file_name_base+'_MSP.txt','w')


# Write the whole array in MASP input format X,Px,Y,Py,Z,Pz,Weight_of_Particle
# This is done in loop - reading each item from Electrons array
# WARNING the momentum is now in p/mc units not SI !!

for i in range(n):                                  
                out.write("%.8e" %(Electrons[i,0]) + \
                " %.8e" %(Electrons[i,1]) + \
                " %.8e" %(Electrons[i,2]) + \
                " %.8e" %(Electrons[i,3]) + \
                " %.8e" %(Electrons[i,4]) + \
                " %.8e" %(Electrons[i,5]) + \
                " %.8e" %(Electrons[i,6]) + "\n")                
                
#Close the files                
out.close()
f.close()
# End of script