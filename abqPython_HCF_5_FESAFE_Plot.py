# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:30:24 2020

@author: tvj
"""

import numpy as np
import matplotlib.pyplot as plt

SultT = 300.E+6
SultC = 1000.E+6
Slife = 150.E+6
filename = 'Job-3-QWosW2-1pt2mm'


"""
LOAD DATA FROM ABAQUS POSTPROCESSING
===============================================================================

FESAFE has 4 columns:
    (1) node ID
    (2) mean stress [MPa]
    (3) stress amplitude [MPa]
    (4) FOS [MPa]
"""
results = np.load(filename + 'Results.odb.npz')
FESAFE = results['FESAFE'].transpose()




"""
CALCULATE GOODMAN RELATIONSHIP AND PLOT SMITH DIAGRAM
===============================================================================
"""
def GoodmanSamp(SultC,SultT,Slife,Smean):
    if Smean >= -1*SultC and Smean < Slife-SultC:
        Samp = SultC * (1 - abs(Smean)/SultC)
    elif Smean >= Slife-SultC and Smean < 0:
        Samp = Slife
    elif Smean >= 0 and Smean <= SultT:
        Samp = Slife * (1 - Smean/SultT)
    else:
        Samp = 0  
    return Samp


# Get the smallest, overall FOS
Goodman_FOS = FESAFE[:,3].min()

# Initialize data points for Smith diagram
Goodman_Smean = np.array([-1*SultC, -1*SultC+Slife, 0, SultT])
Goodman_Samp = np.zeros(len(Goodman_Smean))

# Initialize data points for Smith diagram w/ FOS
Goodman_Smean_FOS = Goodman_Smean / Goodman_FOS
Goodman_Samp_FOS = np.zeros(len(Goodman_Smean_FOS))

for item in range(len(Goodman_Smean)):
    Goodman_Samp[item] = GoodmanSamp(SultC,
                                     SultT,
                                     Slife,
                                     Goodman_Smean[item])
    
    Goodman_Samp_FOS[item] = GoodmanSamp(SultC/Goodman_FOS,
                                         SultT/Goodman_FOS,
                                         Slife/Goodman_FOS,
                                         Goodman_Smean_FOS[item])

Goodman_Smax = Goodman_Smean + Goodman_Samp
Goodman_Smin = Goodman_Smean - Goodman_Samp

Goodman_Smax_FOS = Goodman_Smean_FOS + Goodman_Samp_FOS
Goodman_Smin_FOS = Goodman_Smean_FOS - Goodman_Samp_FOS




"""
PLOT SMITH DIAGRAM
===============================================================================
"""
fig1 = plt.figure(figsize=(6,6), dpi=200)
ax1 = fig1.add_axes([0,0,1,1])
ax1.plot(Goodman_Smean/1.E6,
        Goodman_Smax/1.E6,
        color='blue',
        label="Goodman Max")
ax1.plot(Goodman_Smean/1.E6,
        Goodman_Smin/1.E6,
        color='brown',
        label="Goodman Min")
ax1.plot(Goodman_Smean_FOS/1.E6,
        Goodman_Smax_FOS/1.E6,
        color='blue',
        ls='--',
        label=("{:.2f}".format(Goodman_FOS)+' FOS'))
ax1.plot(Goodman_Smean_FOS/1.E6,
        Goodman_Smin_FOS/1.E6,
        color='brown',
        ls='--',
        label=("(at node "+str(FESAFE[np.argmin(FESAFE[:,3])][0].astype(int))+")"))
ax1.set_xlabel('Mean Stress [MPa]')
ax1.set_title(filename)


# OVERLAY THE FESAFE results (using error bar plotting)

ax1.errorbar(FESAFE[:,1], 
             FESAFE[:,1], 
             yerr = FESAFE[:,2],
             fmt='.',
             markeredgecolor=('k'),
             color='red',
             label="Critical Plane")


plt.legend(loc='lower right')
plt.grid(which='major', axis='both')


fig1.savefig(filename+'_FESAFE_Smith.png', bbox_inches = "tight")



"""
PLOT GOODMAN DIAGRAM
===============================================================================
"""
fig2 = plt.figure(figsize=(6,6), dpi=200)
ax2 = fig2.add_axes([0,0,1,1])

ax2.plot(Goodman_Smean/1.E6,
        Goodman_Samp/1.E6,
        color='blue',
        label="Goodman")
ax2.plot(Goodman_Smean_FOS/1.E6,
        Goodman_Samp_FOS/1.E6,
        color='blue',
        ls='--',
        label=("{:.2f}".format(Goodman_FOS)+' FOS'))
ax2.set_xlabel('Mean Stress [MPa]')
ax2.set_ylabel('Stress Amplitude [MPa]')
ax2.set_title(filename)


# OVERLAY THE FESAFE results

ax2.plot(FESAFE[:,1], 
         FESAFE[:,2], 
         color='red',
         marker='.',
         lw=0,
         label="Critical Plane")

plt.legend(loc='upper right')
plt.grid(which='major', axis='both')

fig2.savefig(filename+'_FESAFE_Goodman.png', bbox_inches = "tight")
