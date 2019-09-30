### psmplr package ###
#
# Posterior distribution sampling framework created to take in an INLA model
# and output a sample to mimic all configurations of the model
#
# Main file containing computationally heavy multiplications and sampling calls
#

require(MASS)
require(Matrix)
require(tidyverse)

# parallel multiplications
selectedMeans <- function(lst_of_means, Amat){
        return(purrr::map(lst_of_means, function(xx) extractEffectMeans(xx, Amat)))
}

selectedCovMat <- function(lst_of_cov, Amat){
        return(purrr::map(lst_of_cov, function(xx) extractEffectCovMat(xx, Amat)))
}

# sampling
psmplr <- function(inla_model, effect_name, n = 1, constraint_point = 2){
        # Part 1 - make the selection matrix
        Amat <- makeAMat(inla_model, effect_name, constraint_point)

        # Part 2 - arguments of MVN sampling function
        n <- sampSizes(inla_model, n)
        mu <- extractAllMeans(inla_model) %>% selectedMeans(Amat)
        Sigma <- extractAllCovMat(inla_model) %>% selectedCovMat(Amat)

        # Part 3 - perform sampling (NOTE: NO FORMATTING DONE YET)
        MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
}
