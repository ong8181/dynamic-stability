# R scripts for Ushio et al. "Fluctuating interaction network and time-varying stability of a natural fish community"

# Liscence
This R code is released under the MIT License, see LICENSE.txt.

# R package
Following packages are required:
rEDM version 0.2.4 (included as a zip file in the inst/ folder), tseriesChaos, deSolve, GGally, intergraph, network, igraph, ggplot2, gridExtra, devtools, pforeach, sna

# Software versions
The original results were generated using rEDM v0.2.4 and R3.2.1 on Mac OSX.

# Installation of the packages
All packages except for pforeach can be installed from CRAN.

To install pforeach from github:
> devtools::install_github("hoxo-m/pforeach")

To install rEDM from the inst/ folder:
> install.packages("inst/rEDM.zip", type = "source", repos = "NULL")

or decompress rEDM.zip and copy it to your R library folder.

# Instructions
The main script file is RunAllScripts.R. Sourcing this file will run the analyses and produce the figures for the paper.
