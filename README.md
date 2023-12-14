# Estimating Symmetries
This repository contains the code and data used in the paper "Estimating Maximal Symmetries of Regression Functions via Subgroup Lattices" by Louis Christie and John Aston. 

The files are split into folders for the relevant sections of the paper. 

Eg 6.1 Contains:
  1) "example_6.1.R" relating to Example 6.1 and the continued appendix. This details a finite case where the search group is the finite dihedral group.
  
Eg 6.2 Contains:
  1) "example_6.2.R" relating to Example 6.2 and the continued appendix. This details the case where the symmetries are certain infinite subgroups of the 3D rotation group SO(3). 
  2) "scenario_1_LCE_perms.csv" containg the simulation results for scenario 1
  3) "scenario_2_LCE_perms.csv" containg the simulation results for scenario 2
  4) "scenario_3_LCE_perms.csv" containg the simulation results for scenario 3
  5) "scenario_4_LCE_perms.csv" containg the simulation results for scenario 4
  6) "scenario_2_LCE_varying_K_2.csv" measuring the effects of different sizes of lattice bases
  7) "scenario_4_LCE_varying_K_2.csv" measuring the effects of different sizes of lattice bases

  All csv files are can be read directly into the R code, or can be generated fresh from the code. 
  
Applications Contains:
  1) "Magneteosphere_and_sunspots.R" used to conduct the testing and estimation for both real world data applications in the paper. 
  2) "SWARM_DATA.csv" containing the data from VIRES used to estimate the symmetry of the magnetic field intensity.
