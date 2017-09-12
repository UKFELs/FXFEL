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
c=2.99792458E8
me=9.11e-31


if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SU2Astra <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()




f=tables.open_file(file_name_in,'r')
Electrons=f.root.Particles.read()

# Descale units from p/mc to SI used by Astra

Electrons[:,1]=Electrons[:,1]*(me*c)
Electrons[:,3]=Electrons[:,3]*(me*c)
Electrons[:,5]=Electrons[:,5]*(me*c)

# Assign number of records (macroparticles) as n
n=len(Electrons)

avg_Z=np.average(Electrons[:,4])
avg_Pz=np.average(Electrons[:,5])
avg_Time=avg_Z/c
avg_chrg=np.average(Electrons[:,6])

refx=0
refy=0
refz=avg_Z
refpx=0
refpy=0
refpz=avg_Pz/5.34428595e-28
refchrg=avg_chrg*e_ch*1.E9
# Create output file named output.txt
out=open(file_name_base+'_SU2A.txt','w')

# create reference particle
out.write("%.8e" %(refx) + " %.8e" %(refy) + " %.8e" %(refz) + " %.8e" %(refpx) + " %.8e" %(refpy) + " %.8e" %(refpz) + " %.8e" %(avg_Time/1.E-9) + " %.8e" %(refchrg) + " %.1i" %(1) + " %.1i" %(5) + "\n") 
#out.write("%.8e" %(refx) + " %.8e" %(refy) + "\n") 

# Write the whole array in MASP input format X,Px,Y,Py,Z,Pz,Weight_of_Particle
# This is done in loop - reading each item from Electrons array



for i in range(n):                                  
                out.write("%.8e" %(Electrons[i,0]) + \
                " %.8e" %(Electrons[i,2]) + \
                " %.8e" %(Electrons[i,4]-avg_Z) + \
                " %.8e" %(Electrons[i,1]/5.34428595e-28) + \
                " %.8e" %(Electrons[i,3]/5.34428595e-28) + \
                " %.8e" %((Electrons[i,5]-avg_Pz)/5.34428595e-28) + \
                " %.8e" %((Electrons[i,4]/c-avg_Time)/1.E-9) + \
                " %.8e" %(Electrons[i,6]*e_ch*1.E9) + \
                " %.1i" %(1) + \
                " %.1i" %(5) + "\n")                        
                
#Close the files                
out.close()
f.close()
# End of script