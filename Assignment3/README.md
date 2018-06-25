Assignment-3 - CFD LAB, Parallelization

GROUP D

Procedure to run code:
-> execute make command in terminal
-> use mpirun -n 'num_proc' ./sim
    num_proc is given by iproc*jproc in the .dat file, and is an integer values

Procedure to visualize:
-> Select the last time step .vtk files for all the sub-domains together and apply.
-> select velocity/glyph for individual sub-domains



Analysis: Lid flow in cavity:
-> the speedup is significantly visible as compared to the sequential code.
-> speedup saturates with increasing number of cores
-> vortex forms below the lid, similar to problem 1. Although the vortex size
  is much smaller as the Reynolds number is much lower.
