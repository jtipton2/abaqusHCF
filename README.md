# abaqusHCF
Using Python to perform High Cycle Fatigue (HCF) analysis of ABAQUS|Explicit transient solutions for brittle materials using Goodman's Method with Maximum Normal Stress Multiaxial Fatigue Failure Theory

## Language:
PYTHON scripting for Abaqus/Viewer Release 2020.HF4


## File Run:
C:\temp>abaqus python abaqusScript_Dynamic_v20200928.py
(Console messages are redirected to abaqus.rpy)


## Summary:

### abqPython_HCF_1_PostProc
This script opens an Abaqus/Explicit solution.  It operates on averaged nodal stresses.  At each timestep, it runs through each node and keeps track of the initial (static), dynamic minimum, and dynamic maximum stresses.  These are saved as python binary files. 

### abqPython_HCF_2_Goodman
This script opens the python binary file.  It calculates a modified Goodman failure envelope and plots a Smith diagram.  For each node of the FEA, it calculates the mean stress and stress amplitude using the peak minimum and maximum stresses over time.  This is repeated for maximum, middle, and minimum principal stresses and is plotted on the Smith diagram.  A Factor of Safety is calculated and plotted using the worst peak stress.  The plot is saved and the data are saved as python binary files.

### abqPython_HCF_3_SaveODB
This script opens the original Abaqus ODB results database.  At the final timestep, it saves the mean stress and stress amplitude for each of the 3 principal stresses as nodal field variables.  It also saves the FOS for each node.


## Source Examples:

### Python

Notes regarding bulkDataBlock:
* https://stackoverflow.com/questions/46573959/accelerate-a-slow-loop-in-abaqus-python-code-for-extracting-strain-data-from-od

Explanations of field output options:
* https://abaqus-docs.mit.edu/2017/English/SIMACAEKERRefMap/simaker-c-fieldoutputcpp.htm

Solution on how to average nodal values:
* https://stackoverflow.com/questions/42783385/max-stress-node/
* https://stackoverflow.com/questions/54350924/find-stresses-at-unique-nodal-on-abaqus-with-python-script

Solution on how to take largest, unaveraged nodal value:
* https://python-forum.io/Thread-Get-max-values-based-on-unique-values-in-another-list-python

### High-Cycle Fatigue

Collins, J. A., Chapter 7 \*High-Cycle Fatigue\* in *Failure of Materials in Mechanical Design*, John Wiley and Sons, 1981.

https://www.tec-science.com/material-science/material-testing/fatigue-limit-diagram-according-to-haigh-and-smith-creation/
