# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:30:24 2020

@author: tvj
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

Sult = 300
Slife = Sult/2
filename = 'Job-3-QWcenterW2-1pt2mm'


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

FOS = np.concatenate((nodeCoord, np.zeros((len(nodeCoord),2))), axis=1)
scratch = np.zeros((3,2))
scratch[0,0] = 1
scratch[1,0] = 2
scratch[2,0] = 3


def getFOS(Sn,Su,prinStress):
    Sm = (prinStress[3] + prinStress[2]) / 2
    Sa = (prinStress[3] - prinStress[2]) / 2
    FOS = 1 / (Sa/Sn + abs(Sm)/Su)
    return FOS


for i in range(len(FOS)):
    scratch[0,1] = getFOS(Slife*1.E6,Sult*1.E6,maxPrin[i])
    scratch[1,1] = getFOS(Slife*1.E6,Sult*1.E6,midPrin[i])
    scratch[2,1] = getFOS(Slife*1.E6,Sult*1.E6,minPrin[i])
    FOS[i,4] = scratch[np.argmin(scratch[:,1])][1]
    FOS[i,5] = scratch[np.argmin(scratch[:,1])][0]




"""
CALCULATE GOODMAN RELATIONSHIP AND PLOT SMITH DIAGRAM
===============================================================================
"""
Goodman_Smean = np.array([-1*Sult, 0, Sult])

Goodman_Samp = Slife * (1 - abs(Goodman_Smean)/Sult)
Goodman_Smax = Goodman_Smean + Goodman_Samp
Goodman_Smin = Goodman_Smean - Goodman_Samp

Goodman_FOS = FOS[:,4].min()

Goodman_Smean_FOS = np.array([-1*Sult/Goodman_FOS, 0, Sult/Goodman_FOS])
Goodman_Samp_FOS = Slife * (1/Goodman_FOS - abs(Goodman_Smean_FOS)/Sult)
Goodman_Smax_FOS = Goodman_Smean_FOS + Goodman_Samp_FOS
Goodman_Smin_FOS = Goodman_Smean_FOS - Goodman_Samp_FOS

fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_axes([0,0,1,1])
ax.plot(Goodman_Smean,
        Goodman_Smax,
        color='blue',
        label="Goodman Max")
ax.plot(Goodman_Smean,
        Goodman_Smin,
        color='brown',
        label="Goodman Min")
ax.plot(Goodman_Smean_FOS,
        Goodman_Smax_FOS,
        color='blue',
        ls='--',
        label=("{:.2f}".format(Goodman_FOS)+' FOS'))
ax.plot(Goodman_Smean_FOS,
        Goodman_Smin_FOS,
        color='brown',
        ls='--',
        label=("(at node "+str(FOS[np.argmin(FOS[:,4])][0].astype(int))+")"))
ax.set_xlabel('Mean Stress [MPa]')
ax.set_title(filename)




"""
OVERLAY THE 3 PRINCIPAL STRESSES (using error bar plotting)
===============================================================================
"""

minPrinAmp = (minPrin[:,3] - minPrin[:,2])/2
minPrinMean = (minPrin[:,3] + minPrin[:,2])/2

ax.errorbar(minPrinMean/1.E6, 
             minPrinMean/1.E6, 
             yerr = minPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='green',
             label="max(S3)")


midPrinAmp = (midPrin[:,3] - midPrin[:,2])/2
midPrinMean = (midPrin[:,3] + midPrin[:,2])/2

ax.errorbar(midPrinMean/1.E6, 
             midPrinMean/1.E6, 
             yerr = midPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='orange',
             label="max(S2)")


maxPrinAmp = (maxPrin[:,3] - maxPrin[:,2])/2
maxPrinMean = (maxPrin[:,3] + maxPrin[:,2])/2

ax.errorbar(maxPrinMean/1.E6, 
             maxPrinMean/1.E6, 
             yerr = maxPrinAmp/1.E6,
             fmt='.',
             markeredgecolor=('k'),
             color='red',
             label="max(S1)")

plt.legend(loc='lower right')
plt.grid(which='major', axis='both')

fig.savefig(filename+'_Goodman.png', bbox_inches = "tight")




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
