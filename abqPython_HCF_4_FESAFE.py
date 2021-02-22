#==============================================================================
# IMPORT NECESSARY MODULES
#==============================================================================
# C:\temp>abaqus viewer -noGUI  (this requires a CAE license)
# C:\temp>abaqus python         (this does not require a license)

import numpy as np
import time
from odbAccess import openOdb
from abaqusConstants import *
from multiprocessing import Pool
from contextlib import closing
import os
import sys




#==============================================================================
# OPEN THE DATABASE
#==============================================================================

filename = 'Job-3-QWosW2-1pt2mmResults.odb'

odb = openOdb(filename,readOnly=True)
Frame = odb.steps['fe-safe_01'].frames

# Get Instance
allInstances = (odb.rootAssembly.instances.keys())
odbInstance = odb.rootAssembly.instances[allInstances[-1]]




#==============================================================================
# GET CRITICAL PLANE NODAL STRESS AMPLITUDES [MPa]
#==============================================================================
Field = Frame[-1].fieldOutputs['Cycle-SAmp'].getSubset(position = ELEMENT_NODAL)

NumValues = len(Field.values)

# Create vector of element nodes and stresses
# (for some reason, ababqus breaks into blocks of data
#  need to join the various data blocks into an array)
Values = Field.bulkDataBlocks[0].data
NodeLabels = Field.bulkDataBlocks[0].nodeLabels
for i in range(len(Field.bulkDataBlocks)-1):
    Values = np.vstack((Values,Field.bulkDataBlocks[i+1].data))
    NodeLabels = np.hstack((NodeLabels,Field.bulkDataBlocks[i+1].nodeLabels))

# Nodes are shared across multiple elements.  Get unique node labels.
NodeLabels_unique, unq_idx = np.unique(NodeLabels, return_inverse=True)
NumNodes = len(NodeLabels_unique)

# Calculate nodal averaged stresses at timestep
Values_Averaged=np.zeros((NodeLabels_unique.size,Values.shape[1]))
unq_counts = np.bincount(unq_idx)
for i in xrange(0,Values.shape[1]):
   ValuesTemp = [item[i] for item in Values]
   unq_sum = np.bincount(unq_idx, weights=ValuesTemp)
   Values_Averaged[:,i] = unq_sum / unq_counts

# Save static stress and also initialize dynamic stress vectors
CriticalAmp = Values_Averaged.copy()




#==============================================================================
# GET CRITICAL PLANE NODAL MEAN STRESS [MPa]
#==============================================================================
Field = Frame[-1].fieldOutputs['Cycle-Sm'].getSubset(position = ELEMENT_NODAL)

NumValues = len(Field.values)

# Create vector of element nodes and stresses
# (for some reason, ababqus breaks into blocks of data
#  need to join the various data blocks into an array)
Values = Field.bulkDataBlocks[0].data
NodeLabels = Field.bulkDataBlocks[0].nodeLabels
for i in range(len(Field.bulkDataBlocks)-1):
    Values = np.vstack((Values,Field.bulkDataBlocks[i+1].data))
    NodeLabels = np.hstack((NodeLabels,Field.bulkDataBlocks[i+1].nodeLabels))

# Nodes are shared across multiple elements.  Get unique node labels.
NodeLabels_unique, unq_idx = np.unique(NodeLabels, return_inverse=True)
NumNodes = len(NodeLabels_unique)

# Calculate nodal averaged stresses at timestep
Values_Averaged=np.zeros((NodeLabels_unique.size,Values.shape[1]))
unq_counts = np.bincount(unq_idx)
for i in xrange(0,Values.shape[1]):
   ValuesTemp = [item[i] for item in Values]
   unq_sum = np.bincount(unq_idx, weights=ValuesTemp)
   Values_Averaged[:,i] = unq_sum / unq_counts

# Save static stress and also initialize dynamic stress vectors
CriticalMean = Values_Averaged.copy()




#==============================================================================
# GET CRITICAL PLANE FACTORS OF SAFETY
#==============================================================================
Field = Frame[-1].fieldOutputs['FOS@Life=Infinite'].getSubset(position = ELEMENT_NODAL)

NumValues = len(Field.values)

# Create vector of element nodes and stresses
# (for some reason, ababqus breaks into blocks of data
#  need to join the various data blocks into an array)
Values = Field.bulkDataBlocks[0].data
NodeLabels = Field.bulkDataBlocks[0].nodeLabels
for i in range(len(Field.bulkDataBlocks)-1):
    Values = np.vstack((Values,Field.bulkDataBlocks[i+1].data))
    NodeLabels = np.hstack((NodeLabels,Field.bulkDataBlocks[i+1].nodeLabels))

# Nodes are shared across multiple elements.  Get unique node labels.
NodeLabels_unique, unq_idx = np.unique(NodeLabels, return_inverse=True)
NumNodes = len(NodeLabels_unique)

# Calculate nodal averaged stresses at timestep
Values_Averaged=np.zeros((NodeLabels_unique.size,Values.shape[1]))
unq_counts = np.bincount(unq_idx)
for i in xrange(0,Values.shape[1]):
   ValuesTemp = [item[i] for item in Values]
   unq_sum = np.bincount(unq_idx, weights=ValuesTemp)
   Values_Averaged[:,i] = unq_sum / unq_counts

# Save static stress and also initialize dynamic stress vectors
CriticalFOS = Values_Averaged.copy()




#==============================================================================
# CLOSE AND SAVE
#==============================================================================
odb.close()

np.savez_compressed(filename, FESAFE=np.vstack((NodeLabels_unique,CriticalMean.flatten(),CriticalAmp.flatten(),CriticalFOS.flatten())))
