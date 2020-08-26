# fdaPDE with automatic optimization methods for GCV function

This repository contains the development version of fdaPDE package. In particular it contains the possibility of using iterative optimization methods (Newton's methods) to find the best smoothing parameter from the minimization of the GCV. 

New features wrt CRAN: smooth regression for manifold and volumetric domains, also with areal data. Smooth fPCA over 2D, 2.5D and 3D domains, also with areal data.

Bugs corrected in the areal data smoothing.

smooth.FEM.basis, smooth.PDE.FEM.basis, smooth.FEM.PDE.sv.basis are deprecated, smooth.FEM has to be used in all cases.

Image.FEM has been restored. Bugs in fPCA, boundary conditions and space-varying regression have been fixed. Issues of point location in 2.5D have been fixed.
Compiled in both Win RStudio and Ubuntu 18.04 using g++ compiler. If using a Linux machine, it is advisable to install rgl, geometry, plot3D and plot3Drgl before fdaPDE. If using Windows, it is advisable to install Rtools and then rgl, plot3D, plot3Drgl, geometry and RcppEigen libraries. 

## Subfolder structure:
/src contains all C++ code and a special file named Makevars necessary to build and install the R package,
/R contains the R functions that wrap the C++ calls,
/data contains the data to run the tests in /tests.
/tests contains the test to be run and the file test_simulations.R, which contains the same test but organized in order to be repeated many times, as used in the report to find dthe boxplopts shown.

## Installation:
Two different methods are proposed in order to install the package in the R environment.  
Download the `.zip` file from the repository, unzip it, and for the installation choose one of the two following methods:  

- R console:
        ```install.packages("/path/to/PACS_MerettPoiatti-master", type='source', repos=NULL)```

- From the Terminal: 
        ```$ R CMD build <path to folder>}```     
        ```$ R CMD INSTALL <path name of the R library tree>}```

see the installation section in the report for more information.
## Remarks:

1) the shift of indexes from R to C++ is done within the R functions smooth.FEM and FPCA.FEM Do not use C++ scripts directly on the R mesh objects, unless you take care of shifing indexes by yourself.
