# psmplr
Source code for psmplr package developed at St. Michael's hospital

This package was inspired by the need to sample the posterior distribution of a random effect of a mixed effect model. When attemping to get predictions from an INLA model we found that it did not store sufficient information about the posterior distribution of the random effects even when predict = TRUE. However, the posterior distribution of the full model is stored and with some considerations it is fairly straight forward to extract the posterior distribution of an effect in order to generate new samples. The thought process and theoretical considerations can be seen on Alin Morariu's website (https://a-morariu.github.io/Spatial_Stats.html).


# How to install

Copy and paste the following: 

library(devtools)
install_github("A-Morariu/psmplr")
