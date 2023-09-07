# biomass yield potential on marginal land
## This repo has all the data and scripts needed to reproduce the results in this paper: _Biomass yield potential on U.S. marginal land and its contribution to reach net-zero emission_ (doi: TBD). The folders are arranged as follows,
- **calibrations**: This includes the optimization, single-site calibration plot and multi-site calibration plot for the three crops. This requires installation of the biocro R package.
- **models**:
  - The BioCro R version used in this paper. Note that it's copied from the biocro-dev repo, and may not be up-to-date with the master branch of biocro-dev.
  - The BioCro Fortran version was used to run the regional simulations. It is currently hosted on [my gitlab repository](https://gitlab.com/yufeng87/biocro_offline_fortran).
- **yield_maps**: This contains some Matlab scripts and related data to plot the yield maps. Running **multi_crop_plot_overlay_yield.m** should reproduce the yield maps shown in the paper.
- **other_plots**: This has one R script to reproduce the biomass bar plot shown in the paper. It also prints out some values to the console, which were used in the paper.
- **marginal_land**: This contains a netcdf file for the marginal land use fractions at CWRF's resolution (30*30 km<sup>2</sup>). 
