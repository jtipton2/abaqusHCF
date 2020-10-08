# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 09:31:17 2020

@author: tvj
"""

import numpy as np
from odbAccess import openOdb
from abaqusConstants import *

filename = 'Job-3-QWcenterW2-1pt2mm'


"""
LOAD DATA
===============================================================================
"""
results = np.load(filename + '.py.npz')
maxPrin = results['maxPrin'].transpose()
minPrin = results['minPrin'].transpose()
midPrin = results['midPrin'].transpose()
nodeCoord = results['nodeCoord']
FOS = results['FOS']




"""
LOAD ODB
===============================================================================
"""

odb = openOdb(filename+'.odb',readOnly=False)


# Get Instance
allInstances = (odb.rootAssembly.instances.keys())
odbInstance = odb.rootAssembly.instances[allInstances[-1]]




"""
FORMAT AND SAVE DATA TO ODB
===============================================================================
"""


# Maximum Principal Stress
SMaxPrinLabels = np.ascontiguousarray(maxPrin[:,0], dtype=np.int)
SMeanMaxPrinValues = np.ascontiguousarray(np.reshape(maxPrin[:,1],(-1,1)), dtype=np.float32)
SAmpMaxPrinValues = np.ascontiguousarray(np.reshape(maxPrin[:,2],(-1,1)), dtype=np.float32)   

newFieldOutputMean = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SMeanMaxPrin', description = 'Mean maximum principal stress', type = SCALAR)
newFieldOutputMean.addData(position=NODAL, instance = odbInstance, labels = SMaxPrinLabels, data = SMeanMaxPrinValues)

newFieldOutputAmp = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SAmpMaxPrin', description = 'Largest maximum principal stress amplitude', type = SCALAR)
newFieldOutputAmp.addData(position=NODAL, instance = odbInstance, labels = SMaxPrinLabels, data = SAmpMaxPrinValues)


# Middle Principal Stress
SMidPrinLabels = np.ascontiguousarray(midPrin[:,0], dtype=np.int)
SMeanMidPrinValues = np.ascontiguousarray(np.reshape(midPrin[:,1],(-1,1)), dtype=np.float32)
SAmpMidPrinValues = np.ascontiguousarray(np.reshape(midPrin[:,2],(-1,1)), dtype=np.float32)   

newFieldOutputMean = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SMeanMidPrin', description = 'Mean middle principal stress', type = SCALAR)
newFieldOutputMean.addData(position=NODAL, instance = odbInstance, labels = SMidPrinLabels, data = SMeanMidPrinValues)

newFieldOutputAmp = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SAmpMidPrin', description = 'Largest middle principal stress amplitude', type = SCALAR)
newFieldOutputAmp.addData(position=NODAL, instance = odbInstance, labels = SMidPrinLabels, data = SAmpMidPrinValues)


# Minimum Principal Stress
SMinPrinLabels = np.ascontiguousarray(minPrin[:,0], dtype=np.int)
SMeanMinPrinValues = np.ascontiguousarray(np.reshape(minPrin[:,1],(-1,1)), dtype=np.float32)
SAmpMinPrinValues = np.ascontiguousarray(np.reshape(minPrin[:,2],(-1,1)), dtype=np.float32)   

newFieldOutputMean = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SMeanMinPrin', description = 'Mean minimum principal stress', type = SCALAR)
newFieldOutputMean.addData(position=NODAL, instance = odbInstance, labels = SMinPrinLabels, data = SMeanMinPrinValues)

newFieldOutputAmp = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'SAmpMinPrin', description = 'Largest minimum principal stress amplitude', type = SCALAR)
newFieldOutputAmp.addData(position=NODAL, instance = odbInstance, labels = SMinPrinLabels, data = SAmpMinPrinValues)


# Goodman Factor of Safety using Brittle Maximum Normal Stress Failure Theory
FOSlabels = np.ascontiguousarray(FOS[:,0], dtype=np.int)
FOSvalues = np.ascontiguousarray(np.reshape(FOS[:,4],(-1,1)), dtype=np.float32)

newFieldOutputFOS = odb.steps['ThermalPulse'].frames[-1].FieldOutput(name = 'FOS', description = 'Maximum normal stress multiaxial fatigue FOS', type = SCALAR)
newFieldOutputFOS.addData(position=NODAL, instance = odbInstance, labels = FOSlabels, data = FOSvalues)




"""
SAVE AND CLOSE
===============================================================================
"""

odb.save()
odb.close()