# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:30:24 2020

@author: tvj
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

SultT = 300.E+6
SultC = 1000.E+6
Slife = 150.E+6
filename = 'Job-3-QWosW2-5mm'


"""
LOAD DATA FROM ABAQUS POSTPROCESSING
===============================================================================

maxPrin, midPrin, and minPrin each have 4 columns:
    (1) node ID
    (2) maximum nodal static stress [Pa]
    (3) minimum nodal dynamic stress [Pa]
    (4) maximum nodal dynamic stress [Pa]
    
nodeCoord has 4 columns:
    (1) node ID
    (2) x coordinate [m]
    (3) y coordinate [m]
    (4) z coordinate [m]
"""
results = np.load(filename + '.odb.npz')
maxPrin = results['maxPrin'].transpose()
minPrin = results['minPrin'].transpose()
midPrin = results['midPrin'].transpose()
nodeCoord = results['nodeCoord']
# Sort nodeCoord on nodal values
nodeCoord = nodeCoord[nodeCoord[:,0].argsort()]




"""
PLOT A POINT CLOUD FOR VISUALIZATION
===============================================================================

plotData = np.concatenate((nodeCoord, maxPrin), axis=1)
# The argsort below sorts the Z coordinate values from min to max
# When the point cloud is plotted, it will layer correctly such that
# the viewer sees only the external surface points
plotData = plotData[plotData[:,3].argsort()[::-1]]

fig = plt.figure(figsize=(10,6), dpi=200)
ax = fig.gca(projection='3d')
p = ax.scatter(plotData[:,1]*100, 
               plotData[:,3]*100, 
               plotData[:,2]*100, 
               c = plotData[:,5]/1.E6,
               depthshade=False,
               cmap=plt.get_cmap("coolwarm"),
               marker=',',
               s = 0.1)

# Make axes limits 
xyzlim = np.array([ax.get_xlim3d(),ax.get_ylim3d(),ax.get_zlim3d()]).T
XYZlim = [min(xyzlim[0]),max(xyzlim[1])]
#XYZlim = [-5,15]
ax.set_xlim3d(XYZlim)
ax.set_ylim3d(XYZlim)
ax.set_zlim3d(XYZlim)

ax.set_xlabel('X (cm)')
ax.set_ylabel('Z (cm)')
ax.set_zlabel('Y (cm)')

#print(ax.azim)
ax.view_init(azim=-120)

#legend
cbar = fig.colorbar(p)
cbar.set_label('Static Max Principal (S1) Stress [MPa]', rotation=270)

plt.show()
"""



"""
CALCULATE FACTORS OF SAFETY
===============================================================================
"""

# Initialize FOS matrix:
# This takes the nodal coordinate matrix and adds 2 empty columns to record
# the minimum FOS at each node along with which principal stress contains
# the minimum.
FOS = np.concatenate((nodeCoord, np.zeros((len(nodeCoord),2))), axis=1)

# Initialze scratch matrix:
# Column [0] is the principal stress direction
# Column [1] will hold the FOS
scratch = np.zeros((3,2))
scratch[0,0] = 1
scratch[1,0] = 2
scratch[2,0] = 3


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


def getFOS(SultC,SultT,Slife,prinStress):
    Smean = (prinStress[3] + prinStress[2]) / 2
    Samp = (prinStress[3] - prinStress[2]) / 2
    
    if Smean >= -1*SultC and Smean < Slife-SultC:
        FOS = 1 / (Samp/SultC + abs(Smean)/SultC)
    elif Smean >= Slife-SultC and Smean < 0:
        FOS = Slife / Samp
    elif Smean >= 0 and Smean <= SultT:
        FOS = 1 / (Samp/Slife + abs(Smean)/SultT)
    else:
        FOS = 0  
    
    return FOS


for i in range(len(FOS)):
    scratch[0,1] = getFOS(SultC,SultT,Slife,maxPrin[i])
    scratch[1,1] = getFOS(SultC,SultT,Slife,midPrin[i])
    scratch[2,1] = getFOS(SultC,SultT,Slife,minPrin[i])
    # Get the smallest FOS
    FOS[i,4] = scratch[np.argmin(scratch[:,1])][1]
    # Get the associated principal stress direction
    FOS[i,5] = scratch[np.argmin(scratch[:,1])][0]



"""
CALCULATE GOODMAN RELATIONSHIP AND PLOT SMITH DIAGRAM
===============================================================================
"""
# Get the smallest, overall FOS
Goodman_FOS = FOS[:,4].min()

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
        label=("(at node "+str(FOS[np.argmin(FOS[:,4])][0].astype(int))+")"))
ax1.set_xlabel('Mean Stress [MPa]')
ax1.set_title(filename)


# OVERLAY THE 3 PRINCIPAL STRESSES (using error bar plotting)

minPrinAmp = (minPrin[:,3] - minPrin[:,2])/2
minPrinMean = (minPrin[:,3] + minPrin[:,2])/2

ax1.errorbar(minPrinMean/1.E6, 
             minPrinMean/1.E6, 
             yerr = minPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='green',
             label="max(S3)")


midPrinAmp = (midPrin[:,3] - midPrin[:,2])/2
midPrinMean = (midPrin[:,3] + midPrin[:,2])/2

ax1.errorbar(midPrinMean/1.E6, 
             midPrinMean/1.E6, 
             yerr = midPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='orange',
             label="max(S2)")


maxPrinAmp = (maxPrin[:,3] - maxPrin[:,2])/2
maxPrinMean = (maxPrin[:,3] + maxPrin[:,2])/2

ax1.errorbar(maxPrinMean/1.E6, 
             maxPrinMean/1.E6, 
             yerr = maxPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='red',
             label="max(S1)")

plt.legend(loc='lower right')
plt.grid(which='major', axis='both')

fig1.savefig(filename+'_Smith.png', bbox_inches = "tight")



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


# PLOT THE 3 PRINCIPAL STRESS AMPLITUDES

ax2.plot(minPrinMean/1.E6, 
         minPrinAmp/1.E6, 
         color='green',
         marker='.',
         lw=0,
         label="max(S3)")

ax2.plot(midPrinMean/1.E6, 
         midPrinAmp/1.E6, 
         color='orange',
         marker='.',
         lw=0,
         label="max(S2)")

ax2.plot(maxPrinMean/1.E6, 
         maxPrinAmp/1.E6, 
         color='red',
         marker='.',
         lw=0,
         label="max(S1)")

plt.legend(loc='upper right')
plt.grid(which='major', axis='both')

fig2.savefig(filename+'_Goodman.png', bbox_inches = "tight")




"""
SAVE DATA
===============================================================================
"""

np.savez_compressed(filename+'.py', 
                    maxPrin=np.vstack((maxPrin[:,0],maxPrinMean,maxPrinAmp)), 
                    midPrin=np.vstack((midPrin[:,0],midPrinMean,midPrinAmp)), 
                    minPrin=np.vstack((minPrin[:,0],minPrinMean,minPrinAmp)), 
                    nodeCoord=nodeCoord, 
                    FOS=FOS)
