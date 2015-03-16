# AWoL-MRF
Label-Fusion Method for MR image data (Validation Code - not a production pipeline) 

Manifest:
* RandomWalksMRI_MRF.m: Main file that runs the experiment
* getNBRInfo: Collects stats and returns seed voxels
* getMSTSequencePlus: returns sequence of voxels based on a minimum spanning tree for a patch
* prim.m: Core spanning tree algorithm (not my code)
