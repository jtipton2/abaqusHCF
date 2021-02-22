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


"""
CALCULATE GOODMAN RELATIONSHIP
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


# Initialize data points for Smith diagram
Goodman_Smean = np.array([-1*SultC, -1*SultC+Slife, 0, SultT])
Goodman_Samp = np.zeros(len(Goodman_Smean))


for item in range(len(Goodman_Smean)):
    Goodman_Samp[item] = GoodmanSamp(SultC,
                                     SultT,
                                     Slife,
                                     Goodman_Smean[item])
    

Goodman_Smax = Goodman_Smean + Goodman_Samp
Goodman_Smin = Goodman_Smean - Goodman_Samp



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
ax1.plot([-1*SultC/1.E6, SultT/1.E6],
        [-1*SultC/1.E6, SultT/1.E6],
        color='black',
        ls='--',
        label="Mean Stress")

ax1.set_xlabel('Mean Stress [MPa]')
ax1.set_ylabel('Minimum or Maximum Stress [MPa]')
ax1.set_title('Estimated Smith Diagram for Irradiated Tungsten')

plt.legend(loc='lower right')
plt.grid(which='major', axis='both')

fig1.savefig('Tungsten_Irradiated_Smith.png', bbox_inches = "tight")



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

ax2.set_xlabel('Mean Stress [MPa]')
ax2.set_ylabel('Stress Amplitude [MPa]')
ax2.set_title('Estimated Goodman Diagram for Irradiated Tungsten')

#plt.legend(loc='upper right')
plt.grid(which='major', axis='both')

fig2.savefig('Tungsten_Irradiated_Goodman.png', bbox_inches = "tight")
