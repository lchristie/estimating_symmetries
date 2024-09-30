# Estimating Symmetries
This repository contains the code and data used in the paper "Estimating Maximal Symmetries of Regression Functions via Subgroup Lattices" by Louis Christie and John Aston. 

This repository contains the following files:

**Figure Generation Files**: 
  - "Figure_5a.R" and "Figure_5b.R" that plot the panels relating to example 6.1 in dimension 2 and 4 in the paper;
  - "Figure_6a.R" and "Figure_6b.R" that plot the panels relating to example 6.2's scenarios 1 and 2;
  - "Figure_7a.R" and "Figure_7b.R" that plot the panels relating to example 6.2's scenarios 3 and 4;
  - "Figure_8.R" that plots both a 3D visualisation of the magnetic field readings and a 2D projection of these measurements;
  - "Figure_9.R" that plots both a 3D visualisation of the sunspot occurances and the lattitude vs time of these occurances;
  - "Figure_10a.R", "Figure_10a.R", and "Figure_10a.R" that plot the continuation of example 6.1 into dimension 6, 10, and 20;
  - "Figure_11a.R" and "Figure_11b.R" that plot the frequency that the cofidence region contains the correct largest symmetry in example 6.1 in dimensions 2 and 4;
  - "Figure_12.R" that plots the effects of changing the lattice base in example 6.2;
  - "Figure_13.R" that draws the Q-Q plots for the distributions on the rotation groups used in the magnetosphere application; and
  - "Figure_14.R" that draws the Q-Q plots for the distributions on the rotation groups used in the sunspot application.

**Table Generation Files**:
  - "Table_1.R" that constructs the data in table 1;
  - "Table_2_Sunspots.R" that constructs the data in table 2; and
  - "Table_3_4_and_5.R" that constructs each of these tables.

**Souce Code Files**:
  - "util_funcs.R" that contains many useful functions for running the tests;
  - "test_function.R" that contains the functions that perform both the asymmetric variation tests and the permutation variant;
  - "ico_syms.R" that containst the axes of the icosohedron and functions that sample uniformly from the gorups of rotation around these axes, and uniformly from SO(3). 

**Data Files**:
  - "SWARM_DATA.csv" that contains the readings of the Magnetic Field.

**Application Files**:
  - "magnetosphere.R" contains the code used to estimate the maximal symmetries of the Earth's magnetic field intensity. No table is produced as all symmetries are rejected.

**Simulation Generation Files**:
  - "example_6.1_generate_sims.R" which generates the simulations for example 6.1 in the paper. The dimension to be tested is changed at the start of the file.
  - "example_6.2_generate_sims.R" which generates the simulations for all four scenarios of example 6.2 in the paper.
  - "example_6.1_high_dim_generate_sims.R" which generates the simulations for the high dimensional example in the appendix. The dimension to be tested is changed at the start of the file.
  - "Example_6.1_Helpers.R" which contains the functions used to compute the confidence region and the estimated largest invariant symmetry in example 6.1 and the high dimensional continuation. 

**Simulation Results**:
  - The files in the /sims/ directory contain the results of running the simulation generation files. These are accesed by the figure generation files, rather than re-running the simulations for every figure generation. 


## How to use this repository

If you wish to rerun any simulations, the process is to first run the relevant simulation generation file (with the appropriate choice of dimension) and then to generate the appropriate figure or table with the relevant file. 

If you wish to use this code for other applications, then the "util_funcs.R" and "test_function.R" can be readily applied to any new data and then the confidence region constructed from the results of these tests. 
