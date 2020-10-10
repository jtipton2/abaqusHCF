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
# DEFINE FUNCTIONS
#==============================================================================


def dynamicAveraged(Frame,StressType):
    
    #
    # INITIALIZE DATA FROM STATIC SOLUTION
    #-------------------------------------
    
    timestep = 0
    
    # Get number of nodes
    Field = Frame[timestep].fieldOutputs['S'].getSubset(position = ELEMENT_NODAL).getScalarField(invariant = StressType)
    # NOTE:  to get a specific set use *.getSubset(region = odbInstance.elementSets['SET-ALL'])
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
    StressStatic = Values_Averaged.copy()
    StressDynamicMin = Values_Averaged.copy()
    StressDynamicMax = Values_Averaged.copy()
    
    
    #
    # CYCLE THROUGH EACH TIMESTEP AND UPDATE DYNAMIC STRESS MIN/MAX DATA
    #-------------------------------------------------------------------
    
    for timestep in range(len(Frame)):
       #
       # Get data from a timestamp
       Field = Frame[timestep].fieldOutputs['S'].getSubset(position = ELEMENT_NODAL).getScalarField(invariant = StressType)
       NumValues = len(Field.values)
       Values = Field.bulkDataBlocks[0].data
       NodeLabels = Field.bulkDataBlocks[0].nodeLabels
       for i in range(len(Field.bulkDataBlocks)-1):
           Values = np.vstack((Values,Field.bulkDataBlocks[i+1].data))
           NodeLabels = np.hstack((NodeLabels,Field.bulkDataBlocks[i+1].nodeLabels))
       #
       # Get average stress at each node
       NodeLabels_unique, unq_idx = np.unique(NodeLabels, return_inverse=True)
       Values_Averaged=np.zeros((NodeLabels_unique.size,Values.shape[1]))
       unq_counts = np.bincount(unq_idx)
       for i in xrange(0,Values.shape[1]):
          ValuesTemp = [item[i] for item in Values]
          unq_sum = np.bincount(unq_idx, weights=ValuesTemp)
          Values_Averaged[:,i] = unq_sum / unq_counts
       max_ind=np.unravel_index(np.argmax(Values_Averaged),Values_Averaged.shape)
       #
       # Update maximum and minimum dynamic stress at each node
       for j in range(NumNodes):
          #
          if Values_Averaged[j] > StressDynamicMax[j]:
             StressDynamicMax[j] = Values_Averaged[j]
             #print("changed max at timestep " + str(timestep) + " and j " + str(j))
          #
          if Values_Averaged[j] < StressDynamicMin[j]:
             StressDynamicMin[j] = Values_Averaged[j]
             #print("changed min at timestep " + str(timestep) + " and j " + str(j))
       # end loop
    # end loop
    
    return NodeLabels_unique, StressStatic.flatten(), StressDynamicMin.flatten(), StressDynamicMax.flatten()




def runAnalysis(filename):
    
    #
    # LOAD ABAQUS SOLUTION DATA
    #-------------------------------------------------------------------
    
    odb = openOdb(filename,readOnly=True)
    Frame = odb.steps['ThermalPulse'].frames
    
    # Get Instance
    allInstances = (odb.rootAssembly.instances.keys())
    odbInstance = odb.rootAssembly.instances[allInstances[-1]]
    
    sys.stdout = open(filename + "_" + str(os.getpid()) + ".out", "w")

    
    sys.stdout.write("\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    sys.stdout.write("\nODB FILENAME = " + filename)
    sys.stdout.write("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    #
    # GET NODAL COORDINATES
    #-------------------------------------------------------------------
    nodeList = Frame[0].fieldOutputs['S'].values[0].instance.nodes
    nodeCoord = np.zeros((len(nodeList),4))
    
    for item in range(len(nodeList)):
        nodeCoord[item,0] = nodeList[item].label
        nodeCoord[item,1] = nodeList[item].coordinates[0]
        nodeCoord[item,2] = nodeList[item].coordinates[1]
        nodeCoord[item,3] = nodeList[item].coordinates[2]
    
    #
    # PROCESS RESULTS
    #-------------------------------------------------------------------
    sys.stdout.write("\n\nMAXIMUM PRINCIPAL STRESS")
    sys.stdout.write("\n========================")
    start_time = time.time()
    maxPrin = dynamicAveraged(Frame,MAX_PRINCIPAL)
    sys.stdout.write("\nThe number of unique nodes is " +str(len(maxPrin[0])))
    sys.stdout.write("\nThe max static stress is at NodeLabel "+str(maxPrin[0][np.argmax(maxPrin[1])].astype(int))+ " its value is "+ str(maxPrin[1].max()) +" Pa.")
    StressDynamicAmp = (maxPrin[3] - maxPrin[2]) / 2
    sys.stdout.write("\nThe max stress amplitude is at NodeLabel "+str(maxPrin[0][np.argmax(StressDynamicAmp)].astype(int))+ " its value is "+ str(StressDynamicAmp.max()) +" Pa.")
    sys.stdout.write("\n--- %s seconds ---" % (time.time() - start_time))
    
    
    sys.stdout.write("\n\nMIDDLE PRINCIPAL STRESS")
    sys.stdout.write("\n========================")
    start_time = time.time()
    midPrin = dynamicAveraged(Frame,MID_PRINCIPAL)
    sys.stdout.write("\nThe number of unique nodes is " +str(len(midPrin[0])))
    sys.stdout.write("\nThe max static stress is at NodeLabel "+str(midPrin[0][np.argmax(midPrin[1])].astype(int))+ " its value is "+ str(midPrin[1].max()) +" Pa.")
    StressDynamicAmp = (midPrin[3] - midPrin[2]) / 2
    sys.stdout.write("\nThe max stress amplitude is at NodeLabel "+str(midPrin[0][np.argmax(StressDynamicAmp)].astype(int))+ " its value is "+ str(StressDynamicAmp.max()) +" Pa.")
    sys.stdout.write("\n--- %s seconds ---" % (time.time() - start_time))
    
    
    sys.stdout.write("\n\nMINIMUM PRINCIPAL STRESS")
    sys.stdout.write("\n========================")
    start_time = time.time()
    minPrin = dynamicAveraged(Frame,MIN_PRINCIPAL)
    sys.stdout.write("\nThe number of unique nodes is " +str(len(minPrin[0])))
    sys.stdout.write("\nThe max static stress is at NodeLabel "+str(minPrin[0][np.argmax(minPrin[1])].astype(int))+ " its value is "+ str(minPrin[1].max()) +" Pa.")
    StressDynamicAmp = (minPrin[3] - minPrin[2]) / 2
    sys.stdout.write("\nThe max stress amplitude is at NodeLabel "+str(minPrin[0][np.argmax(StressDynamicAmp)].astype(int))+ " its value is "+ str(StressDynamicAmp.max()) +" Pa.")
    sys.stdout.write("\n--- %s seconds ---" % (time.time() - start_time))
    
    sys.stdout.close()
   
    
    #
    # SAVE DATA TO NUMPY COMPRESSED BINARY FILE FOR LATER USE
    #-------------------------------------------------------------------
    
    np.savez_compressed(filename, maxPrin=maxPrin, midPrin=midPrin, minPrin=minPrin, nodeCoord=nodeCoord)
    
    # results = np.load('filename.npx')
    # results['maxPrin']
    
    
    #
    # TERMINATE PROGRAM
    #-------------------------------------------------------------------
    
    #odb.save()
    odb.close()




#==============================================================================
# RUN THE PROGRAM
#==============================================================================

# Consider running in parallel for batch jobs:
# <https://stackabuse.com/parallel-processing-in-python/>
# <https://stackoverflow.com/a/7207336>
# <https://stackoverflow.com/a/25968716/1426569>


filenames = ['Job-3-QWosW2-5mm.odb',
             'Job-3-QWosW2-3pt5mm.odb']
#             'Job-3-QWosW2-2pt45mm.odb',
#             'Job-3-QWosW2-1pt715mm.odb',
#             'Job-3-QWosW2-1pt2mm.odb',
#             'Job-3-QWcenterW2-1pt2mm.odb']
#             'Job-3-QWosW2-0pt84mm.odb']

if __name__ == '__main__':
    agents = 6
    with closing(Pool(processes=agents)) as pool:
        pool.map(runAnalysis, filenames)
        pool.terminate()

