# bacterial_cellwall_division

To install a Fortran compiler in Ubuntu, run the following in a terminal:
sudo apt-get install gfortran

To compile REMODELER 2, run the following in a terminal (which takes ~ 1 min):

gfortran -fopenmp -fno-range-check Vars.f90 Mods.f90 Remodeler2.f90 -o wall

To run simulation:

./wall >log.txt


Note: the software has been tested on Linux operating systems including Ubuntu 16.04.4 LTS

The demo simulation tests the effect of PG remodeling in the presence of a constrictive force. Demo output: visual.psf and visual001.dcd

The data can be visualized using VMD, which can be downloaded from: https://www.ks.uiuc.edu/Research/vmd/

in a terminal, run the following to visualize the output:
vmd visual.psf and visual001.dcd

