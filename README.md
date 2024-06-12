# water_shuffle
Permute waters in an MD trajectory to minimise diffusion

## Motivation
For the analysis of hydration patters round solute molecules, it is sometimes
useful to permute the indices of the water molecules in each frame of an MD
trajectory so that each then appears to remain within a restricted volume.

## Method
The Python code here implements this using the linear sum assigment method in a iterative manner:

1. The coordinates of each water oxygen atom in each trajectory frame are selected.
2. The coordinates in the first frame are copied to form the initial guess for the mean coordinates.
3. For each frame, the distance matrix of the water oxygen coordinates from those of the mean is calculated. The linear sum assignment method (as implemented in scipy) is used to find a
permutation vector that is then used to reorder the water oxygen atoms i(and associated hydrogens) to minimise the distances of each from the corresponding mean. 
5. Once each frame has been processed, the mean coordinates of each water molecule (over all trajectory frames) is updated.
7. Until a set number of iterations is exceeded, the process is repeated from step 3.
8. The permuted trajectory is saved to a new file.

## Installation

### Requirements
* Python (version 3.9 or later)
* git
* pip

The software has been tested on Linux (Ubuntu) and Mac (Ventura) platforms, with Python versions 3.9, 3.10 and 3.11.

### Procedure
In a suitable directory:

```
git clone https://github.com/CharlieLaughton/water_shuffle.git
cd water_shuffle
pip install .
```
Installation will typically take less than one minute. You should then find the command `water-shuffle` in your path:
```
which water-shuffle`
```

### Testing
The `/test` directory contains a small trajectory and topology file for testing (alanine dipeptide in a box of water).

To run:
```
water-shuffle --trajin test.nc --topology test.prmtop --water_name HOH --n_cycles 5 --trajout permuted.nc
```
The test job should take less than 30 seconds to run. If you compare the motion of a single water molecule (e.g. residue 4) in the original and permuted versions of the trajectory file (e.g. using VMD or similar) the removal of diffusive motion will be apparent.
