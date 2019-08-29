# psmplr package
#
# Posterior distribution sampling framework created to take in an INLA model
# and output a sample to mimic all configurations of the model
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

require(MASS)
require(Matrix)
require(tidyverse)

reID <- function(inla_model, re_name){ ### FIX TO WHICH INSTEAD OF GREP
        ### return a vector of the indicies (row and column) of the sub-matrix
        ### corresponding to the named random effect (ONLY works for 1 effect)
        name <- grep(re_name,inla_model$misc$configs$contents$tag)
        indicies <- seq(inla_model$misc$configs$contents$start[name],
                len = inla_model$misc$configs$contents$length[name])

        return(indicies)
}

createTransform <- function(re_index, inla_model, constraint_point){
        ### return the A matrix which selects the sub-matrix corresponding to
        ### the named random effect while placing a constraint at a single
        ### point in the walk
        # Take the sub-matrix of Q that we want
        # number of columns - length of indices (for the random effect of interest)
        # number of rows - to match dimensions of Q (for multiplication)

        dimension <- nrow(inla_model$misc$configs$config[[1]]$Q)

        ### Top
        top_block <- Matrix::Matrix(0,
                nrow = (min(re_index)-1),
                ncol = length(re_index))

        ### Bottom
        bot_block <- Matrix::Matrix(0,
                nrow = (dimension - max(re_index)),
                ncol = length(re_index))

        ### Middle
        if(is.null(constraint_point)){
                mid_block <- Matrix::Diagonal(n = max(re_index) - min(re_index) + 1)
        } else{
                ### Middle
                differencing_mat <- Matrix(0,
                        nrow = max(re_index) - min(re_index) + 1,
                        ncol = max(re_index) - min(re_index) + 1)
                differencing_mat[,as.numeric(constraint_point)] <- 1

                # get +/- 1 entries on each row except in at entry
                # [fixed_point, constraint_point] where we get a row of zeros
                mid_block <- Matrix::Diagonal(n = max(re_index) - min(re_index) + 1) - differencing_mat
        }
        return(rbind(top_block, mid_block, bot_block))
}

##### Helper functions #####

# selection matrix
makeAMat <- function(inla_model, effect_name, contraint_point = 2){
        return(inla_model %>%
                reID(effect_name) %>%
                createTransform(inla_model, contraint_point))
}

# sample size manipulation
s.weights <- function(inla_model){
        weights <- c()
        for (i in 1:inla_model$misc$configs$nconfig) {
                weights[i] <- inla_model$misc$configs$config[[i]]$log.posterior
        }
        return(weights/sum(weights))
}

sampSizes <- function(inla_model, n = 1){
        return(as.vector(ceiling(s.weights(inla_model)*n) ) )
}

# extract the necessary components for sampling from INLA model
extractAllMeans <- function(inla_model){
        return(purrr::map(inla_model$misc$configs$config, function(xx) xx$mean))
}

extractAllCovMat <- function(inla_model){
        return(purrr::map(inla_model$misc$configs$config, function(xx) xx$Q))
}

# keep heavy multiplication steps seperate
extractEffectMeans <- function(mu, Amat){
        # takes in one mean vector and selects the bits we need based on the AMat
        return(as.vector(Matrix::crossprod(Amat,as.matrix(mu) ) ) )
}

extractEffectCovMat <- function(sigma, Amat){
        # takes in one precision matrix and selects the submatrix we need based on the AMat
        return(new("dsCMatrix",
                x = sigma@x,
                i = sigma@i,
                p = sigma@p,
                Dim = sigma@Dim) %>%
                Matrix::Cholesky(LDL = FALSE, perm = FALSE) %>%
                Matrix::solve(Amat) %>%
                Matrix::crossprod() )
}

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
