# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 10:47:58 2017

@author: ptracz
"""

import numpy as np
import tables
import sys
import matplotlib.pyplot as plt
c=3.0e+8                    # Speed of light
m=9.11e-31                  # mass of electron
e_ch=1.602e-19

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: Emittance <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()
 
f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()

Pi=np.pi                    # Pi number taken from 'numpy' as more precise than just 3.1415
c=3.0e+8                    # Speed of light
m=9.11e-31                  # mass of electron
e_ch=1.602e-19              # charge of electron
# Sort particles along Z-axis
Particles=Particles[Particles[:,4].argsort()]


def Weighted_Stats(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2.0, weights=weights)  # Fast and numerically precise
    return (np.sqrt(variance))

print 'std.Dev wght X =',Weighted_Stats(Particles[:,0], Particles[:,6])
print 'std.Dev wght Px =',Weighted_Stats(Particles[:,1], Particles[:,6])
print 'std.Dev wght Y =',Weighted_Stats(Particles[:,2], Particles[:,6])
print 'std.Dev wght Py =',Weighted_Stats(Particles[:,3], Particles[:,6])
print 'std.Dev wght Z =',Weighted_Stats(Particles[:,4], Particles[:,6])
print 'std.Dev wght Pz =',Weighted_Stats(Particles[:,5], Particles[:,6])
print 'std.Dev wght NE =',Weighted_Stats(Particles[:,6], Particles[:,6])

# Rescale to SI units
Particles[:,1]=Particles[:,1]*(m*c)
Particles[:,3]=Particles[:,3]*(m*c)
Particles[:,5]=Particles[:,5]*(m*c)

print 'Total charge of the beam = ',np.sum(Particles[:,6])*e_ch

#==============================================================================
# plt.hist(Particles[:,5], bins=40,weights=Particles[:,6])  # arguments are passed to np.histogram
# plt.title("Histogram Z")
# plt.show()
#==============================================================================

def Calculate_Emittance(i):
    z_low=np.min(Particles[:,4])+(i*Step_Z)
    z_high=np.min(Particles[:,4])+((i+1)*Step_Z)
    z_pos=(z_high+z_low)/2.0
    
    m_MOMz=Particles[:,5][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
    
    m_POSx=Particles[:,0][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
    m_MOMx=Particles[:,1][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
    
    m_POSy=Particles[:,2][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
    m_MOMy=Particles[:,3][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]    
    
    x_2=((np.sum(m_POSx*m_POSx))/len(m_POSx))-(np.mean(m_POSx))**2.0
    px_2=((np.sum(m_MOMx*m_MOMx))/len(m_MOMx))-(np.mean(m_MOMx))**2.0
    xpx=np.sum(m_POSx*m_MOMx)/len(m_POSx)-np.sum(m_POSx)*np.sum(m_MOMx)/(len(m_POSx))**2.0
    
    y_2=((np.sum(m_POSy*m_POSy))/len(m_POSy))-(np.mean(m_POSy))**2.0
    py_2=((np.sum(m_MOMy*m_MOMy))/len(m_MOMy))-(np.mean(m_MOMy))**2.0
    ypy=np.sum(m_POSy*m_MOMy)/len(m_POSy)-np.sum(m_POSy)*np.sum(m_MOMy)/(len(m_POSy))**2.0
    
    eps_rms_x=(1.0/(m*c))*np.sqrt((x_2*px_2)-(xpx*xpx))
    eps_rms_y=(1.0/(m*c))*np.sqrt((y_2*py_2)-(ypy*ypy))
    
    p_tot_slc=np.sqrt((m_MOMx[:]**2.0)+(m_MOMy[:]**2.0)+(m_MOMz[:]**2.0))
    gamma_slc=(np.sqrt(1+(p_tot_slc/(m*c))**2.0))
    charge_slc=(np.sum(Particles[:,6][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]))*e_ch
    current_slc=(charge_slc/(z_high-z_low))*c
#    print 'Emit_x = ',eps_rms_x,'Emit_y = ',eps_rms_y,'Z = ',z_pos
    return z_pos,eps_rms_x,eps_rms_y,np.mean(gamma_slc),np.std(gamma_slc),charge_slc,current_slc
    
Num_Slices=100
Step_Z=(np.max(Particles[:,4])-np.min(Particles[:,4]))/Num_Slices


import multiprocessing
pool = multiprocessing.Pool()
if __name__ == '__main__':
    emit = []
    for i in range(0,Num_Slices):
        z_low=np.min(Particles[:,4])+(i*Step_Z)
        z_high=np.min(Particles[:,4])+((i+1)*Step_Z)
        m_MOMz=Particles[:,5][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
        if (len(m_MOMz)>3):
            pool.apply_async(Calculate_Emittance, (i,),callback=emit.append)
    pool.close()
    pool.join()

Num_Slices=len(emit)
full_emit=np.zeros((Num_Slices,3))
full_gamma=np.zeros((Num_Slices,3))
full_charge=np.zeros(Num_Slices)
full_current=np.zeros(Num_Slices)

for i in range(0,Num_Slices):
    full_emit[i,0]=emit[i][0]
    full_emit[i,1]=emit[i][1]
    full_emit[i,2]=emit[i][2]
# The below value is Gamma !!!
    full_gamma[i,0]=emit[i][0]
    full_gamma[i,1]=emit[i][3]
# The below value is Delta_Gamma !!!
    full_gamma[i,2]=emit[i][4]
    full_charge[i]=emit[i][5]
    full_current[i]=emit[i][6]
full_emit=full_emit[full_emit[:,0].argsort()]
full_gamma=full_gamma[full_gamma[:,0].argsort()]  
print 'Total charge of slices = ', np.sum(full_charge) 

from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick
fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, sharey=False, figsize=(10,10))

#plot1=fig.add_subplot(2,2,1)

#ax[0].set_xlabel('Z')
ax[0].set_ylabel('Emittance')
ax[0].yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
#ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
ax[0].plot(full_emit[:,0],full_emit[:,1],label='Emittance X')
ax[0].plot(full_emit[:,0],full_emit[:,2],label='Emittance Y')
ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax[1].set_ylabel(r'$\gamma$')
ax[1].yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
ax[1].plot(full_gamma[:,0],full_gamma[:,1])
ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax[2].set_ylabel(r'$\delta \gamma$')
ax[2].yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
ax[2].plot(full_gamma[:,0],full_gamma[:,2])
ax[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax[3].plot(full_gamma[:,0],full_current)
ax[3].set_xlabel('Z')
ax[3].yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
ax[3].xaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
ax[3].set_ylabel('Current [A]')
ax[3].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax[3].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

plt.savefig(file_name_base+'_analysis.png')
plt.show()






